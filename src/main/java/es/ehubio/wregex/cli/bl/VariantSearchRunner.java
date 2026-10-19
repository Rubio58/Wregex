package es.ehubio.wregex.cli.bl;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.io.Writer;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.zip.GZIPInputStream;

import es.ehubio.db.fasta.Fasta;
import es.ehubio.wregex.Pssm;
import es.ehubio.wregex.Result;
import es.ehubio.wregex.ResultGroup;
import es.ehubio.wregex.Wregex;
import es.ehubio.wregex.cli.dl.HgvsIndex;
import es.ehubio.wregex.cli.dl.MotifCatalog;
import es.ehubio.wregex.cli.dl.MotifDefinition;
import es.ehubio.wregex.cli.dl.MotifInformation;
import es.ehubio.wregex.cli.dl.ProteinEntry;
import es.ehubio.wregex.cli.dl.ProteinIndex;
import es.ehubio.wregex.cli.dl.ProteinRef;
import es.ehubio.wregex.cli.dl.UniProtMap;
import es.ehubio.wregex.cli.dl.VariantRecord;
import es.ehubio.wregex.cli.dl.VariantTsvReader;

public final class VariantSearchRunner {

    // Rutas externas hardcodeadas (temporalmente)
    private static final String HGVS_PATH      = "/home/rubio/Downloads/hgvs4variation.txt.gz";
    private static final String IDMAPPING_PATH = "/home/rubio/Downloads/HUMAN_9606_idmapping_selected.tab.gz";
    private static final String REFSEQ_DIR     = "/home/rubio/Downloads/refseq_human_proteins";

    // TSV embebido en el classpath
    private static final String TSV_RESOURCE = "/clinvar_download_20260924_131330.tsv";

    private static final int FLANK = 20;

    private final String motif;
    private final File   outputDir;

    public VariantSearchRunner(String motif, File outputDir) {
        this.motif = motif;
        this.outputDir = outputDir;
    }

    public List<VariantResult> run() throws Exception {
        Wregex wregex = buildWregex(motif);

        // 1. TSV
        List<VariantRecord> records;
        try (InputStream in = openClasspath(TSV_RESOURCE);
             Reader r = new InputStreamReader(in, StandardCharsets.UTF_8)) {
            records = VariantTsvReader.load(r);
        }

        // 2. VariationIDs únicos
        Set<String> variationIds = new HashSet<>();
        for (VariantRecord r : records) variationIds.add(r.getVariationId());

        // 3. hgvs4variation filtrado + descartes
        File discardedFile = new File(outputDir,
                "discarded_no_protein_expression_" + timestamp() + ".tsv");
        HgvsIndex hgvs;
        try (Reader r = new InputStreamReader(
                    new GZIPInputStream(new FileInputStream(HGVS_PATH)),
                    StandardCharsets.UTF_8);
             Writer discarded = new FileWriter(discardedFile)) {
            discarded.write("symbol\tvariationId\tnucleotideExpression\tproteinExpression\n");
            hgvs = HgvsIndex.load(r, variationIds, discarded);
        }
        System.err.println("[VariantSearchRunner] discarded file: " + discardedFile.getName());

        // 4. proteinIDs únicos (sin versión)
        Set<String> proteinIds = new HashSet<>();
        for (VariantRecord r : records)
            for (ProteinRef ref : hgvs.findByVariationId(r.getVariationId()))
                proteinIds.add(ref.getProteinId());

        // 5. UniProtMap
        UniProtMap uniprotMap;
        try (Reader r = new InputStreamReader(
                new GZIPInputStream(new FileInputStream(IDMAPPING_PATH)),
                StandardCharsets.UTF_8)) {
            uniprotMap = UniProtMap.load(r, proteinIds);
        }

        // 6. FASTA RefSeq filtrado
        ProteinIndex index = ProteinIndex.loadDirectory(new File(REFSEQ_DIR), proteinIds);

        // 7. Enriquecer ProteinRef con uniProt y proteinSize
        hgvs.enrich(uniprotMap, index);

        // 8. Iterar
        List<VariantResult> results = new ArrayList<>();
        for (VariantRecord record : records) {
            List<ProteinRef> refs = hgvs.findByVariationId(record.getVariationId());
            if (refs.isEmpty()) {
                results.add(new VariantResult(record,
                        Collections.<ProteinRef>emptyList(),
                        "SKIPPED", "NO_HGVS", "", "", null, null));
                continue;
            }
            for (ProteinRef ref : refs) {
                VariantResult vr = process(wregex, index, record, ref);
                if (vr != null) results.add(vr);
            }
        }

        System.err.println("[VariantSearchRunner] processed " + results.size() + " rows");
        return results;
    }

    private VariantResult process(Wregex wregex, ProteinIndex index,
                                  VariantRecord record, ProteinRef ref) {
        ProteinEntry protein = index.findByProteinId(ref.getProteinId());
        if (protein == null) return null;

        int pos = ref.getPosition();
        if (pos > protein.length())
            return single(record, ref, "SKIPPED", "OUT_OF_RANGE", "", "", null, null);

        if (protein.aaAt(pos) != ref.getWtAA())
            return single(record, ref, "SKIPPED", "WT_MISMATCH", "", "", null, null);

        int start = Math.max(1, pos - FLANK);
        int end   = Math.min(protein.length(), pos + FLANK);
        String wt = protein.getSequence().substring(start - 1, end);

        int centerOffset = pos - start;
        char[] chars = wt.toCharArray();
        chars[centerOffset] = ref.getMutAA();
        String mut = new String(chars);

        Result wtMatch, mutMatch;
        try {
            wtMatch  = searchBest(wregex, wt, record, ref);
            mutMatch = searchBest(wregex, mut, record, ref);
        } catch (Fasta.InvalidSequenceException e) {
            return single(record, ref, "SKIPPED", "WINDOW_INVALID", wt, mut, null, null);
        }

        return single(record, ref, "RESOLVED", "", wt, mut, wtMatch, mutMatch);
    }

    private static VariantResult single(VariantRecord record, ProteinRef ref,
                                        String status, String reason,
                                        String wt, String mut,
                                        Result wtMatch, Result mutMatch) {
        return new VariantResult(record,
                ref == null ? Collections.<ProteinRef>emptyList() : Collections.singletonList(ref),
                status, reason, wt, mut, wtMatch, mutMatch);
    }

    private Result searchBest(Wregex wregex, String sequence,
                              VariantRecord record, ProteinRef ref)
            throws Fasta.InvalidSequenceException {
        Fasta fasta = new Fasta(record.getName() + "|" + ref, sequence,
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

    private static InputStream openClasspath(String resource) {
        InputStream in = VariantSearchRunner.class.getResourceAsStream(resource);
        if (in == null)
            throw new IllegalStateException("Missing bundled resource: " + resource);
        return in;
    }

    private static String timestamp() {
        return new java.text.SimpleDateFormat("yyyyMMdd_HHmmss").format(new java.util.Date());
    }
}