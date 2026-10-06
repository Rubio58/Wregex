package es.ehubio.wregex.cli.bl;

import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPInputStream;

import es.ehubio.db.fasta.Fasta;
import es.ehubio.wregex.Pssm;
import es.ehubio.wregex.Result;
import es.ehubio.wregex.ResultGroup;
import es.ehubio.wregex.Wregex;
import es.ehubio.wregex.cli.dl.MotifCatalog;
import es.ehubio.wregex.cli.dl.MotifDefinition;
import es.ehubio.wregex.cli.dl.MotifInformation;
import es.ehubio.wregex.cli.dl.ProteinChange;
import es.ehubio.wregex.cli.dl.ProteinEntry;
import es.ehubio.wregex.cli.dl.ProteinIndex;
import es.ehubio.wregex.cli.dl.VariantRecord;
import es.ehubio.wregex.cli.dl.VariantTsvReader;

public final class VariantSearchRunner {

    private static final String TSV_RESOURCE   = "/clinvar_download_20260924_131330.tsv";
    private static final String FASTA_RESOURCE = "/uniprotkb_AND_model_organism_9606_2026_09_24.fasta.gz";
    private static final int    FLANK          = 10;

    private final String motif;

    public VariantSearchRunner(String motif) {
        this.motif = motif;
    }

    public List<VariantResult> run() throws Exception {
        Wregex wregex = buildWregex(motif);

        List<VariantRecord> records;
        try (InputStream in = open(TSV_RESOURCE);
             Reader r = new InputStreamReader(in, StandardCharsets.UTF_8)) {
            records = VariantTsvReader.load(r);
        }

        ProteinIndex index;
        try (InputStream in = open(FASTA_RESOURCE);
             Reader r = new InputStreamReader(new GZIPInputStream(in), StandardCharsets.UTF_8)) {
            index = ProteinIndex.load(r);
        }

        List<VariantResult> results = new ArrayList<>();
        for (VariantRecord record : records)
            for (ProteinChange change : record.getProteinChanges())
                results.add(process(wregex, index, record, change));

        System.err.println("[VariantSearchRunner] processed " + results.size() + " mutations");
        return results;
    }

    private VariantResult process(Wregex wregex, ProteinIndex index,
                                  VariantRecord record, ProteinChange change) {
        // 1. Resolver proteína
        ProteinEntry protein = null;
        String usedGene = null;
        for (String gene : record.getGeneCandidates()) {
            for (ProteinEntry p : index.findByGene(gene)) {
                if (change.getPosition() <= p.length()
                        && p.aaAt(change.getPosition()) == change.getWtAA()) {
                    protein = p;
                    usedGene = gene;
                    break;
                }
            }
            if (protein != null) break;
        }
        if (protein == null)
            return new VariantResult(record, change, "SKIPPED", "NO_PROTEIN",
                    "", "", "", "", null, null);

        // 2. Extraer ventana WT y mutada
        int pos = change.getPosition();
        int len = protein.length();
        int start = Math.max(1, pos - FLANK);
        int end   = Math.min(len, pos + FLANK);
        String wt = protein.getSequence().substring(start - 1, end);

        int centerOffset = pos - start;
        char[] chars = wt.toCharArray();
        chars[centerOffset] = change.getMutAA();
        String mut = new String(chars);

        // 3. Buscar en WT y en mutada
        Result wtMatch, mutMatch;
        try {
            wtMatch  = searchBest(wregex, wt,  record, change);
            mutMatch = searchBest(wregex, mut, record, change);
        } catch (Fasta.InvalidSequenceException e) {
            return new VariantResult(record, change, "SKIPPED", "WINDOW_INVALID",
                    usedGene, protein.getAccession(), wt, mut, null, null);
        }

        return new VariantResult(record, change, "RESOLVED", "",
                usedGene, protein.getAccession(), wt, mut, wtMatch, mutMatch);
    }

    private Result searchBest(Wregex wregex, String sequence,
                              VariantRecord record, ProteinChange change)
            throws Fasta.InvalidSequenceException {
        Fasta fasta = new Fasta(record.getName() + "|" + change, sequence,
                Fasta.SequenceType.PROTEIN);
        Result best = null;
        for (ResultGroup g : wregex.searchGrouping(fasta))
            for (Result r : g)
                if (best == null || r.getScore() > best.getScore())
                    best = r;
        return best;
    }

    private Wregex buildWregex(String token) throws Exception {
        MotifCatalog catalog = MotifCatalog.load();
        String[] parts = token.split("@", 2);
        String motifName = parts[0];
        String defName = parts.length > 1 ? parts[1] : null;

        MotifInformation motif = catalog.findMotif(motifName);
        if (motif == null)
            throw new IllegalArgumentException("Unknown motif '" + motifName + "'");

        MotifDefinition def = defName == null
                ? motif.getDefinitions().get(0)
                : motif.getDefinitions().stream()
                        .filter(d -> d.getName().equalsIgnoreCase(defName))
                        .findFirst()
                        .orElseThrow(() -> new IllegalArgumentException(
                                "Unknown definition '" + defName + "' for '" + motifName + "'"));

        Pssm pssm = def.getPssm() == null ? null : catalog.loadPssm(def.getPssm());
        return new Wregex(def.getRegex(), pssm);
    }

    private static InputStream open(String resource) {
        InputStream in = VariantSearchRunner.class.getResourceAsStream(resource);
        if (in == null)
            throw new IllegalStateException("Missing bundled resource: " + resource);
        return in;
    }
}