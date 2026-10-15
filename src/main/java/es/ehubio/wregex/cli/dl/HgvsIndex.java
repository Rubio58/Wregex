package es.ehubio.wregex.cli.dl;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Reader;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public final class HgvsIndex {

    // Formato: NP_000654.2:p.Phe213Cys
    private static final Pattern PROTEIN_EXPR =
            Pattern.compile("^(NP_\\d+)\\.(\\d+):p\\.([A-Z][a-z]{2})(\\d+)([A-Z][a-z]{2})$");

    private static final int COL_SYMBOL       = 0;
    private static final int COL_VARIATION_ID = 2;
    private static final int COL_TYPE         = 4;
    private static final int COL_NUCLEOTIDE   = 6;
    private static final int COL_PROTEIN_EXPR = 8;

    private final Map<String, List<ProteinRef>> byVariationId;

    private HgvsIndex(Map<String, List<ProteinRef>> byVariationId) {
        this.byVariationId = byVariationId;
    }

    /**
     * @param reader           lector del hgvs4variation.txt.gz (ya descomprimido)
     * @param wantedVariationIds conjunto de VariationIDs que nos interesan
     * @param discardedOut     writer donde volcar las entradas descartadas
     *                         (puede ser null si no se quiere escribir)
     */
    public static HgvsIndex load(Reader reader,
                                 Set<String> wantedVariationIds,
                                 Writer discardedOut) throws IOException {
        Map<String, List<ProteinRef>> byVariationId = new HashMap<>();
        Set<String> seen = new HashSet<>();
        BufferedReader br = new BufferedReader(reader);
        String line;
        boolean headerSkipped = false;

        while ((line = br.readLine()) != null) {
            if (!headerSkipped) {
                if (line.startsWith("#Symbol")) headerSkipped = true;
                continue;
            }
            if (line.isBlank()) continue;

            String[] c = line.split("\t", -1);
            if (c.length < 9) continue;

            String variationId = c[COL_VARIATION_ID];
            if (!wantedVariationIds.contains(variationId)) continue;
            if (!"coding".equals(c[COL_TYPE])) continue;

            String proteinExpr = c[COL_PROTEIN_EXPR];
            if (proteinExpr == null || proteinExpr.equals("-") || proteinExpr.isEmpty()) {
                writeDiscarded(discardedOut, c);
                continue;
            }

            Matcher m = PROTEIN_EXPR.matcher(proteinExpr);
            if (!m.matches()) {
                writeDiscarded(discardedOut, c);
                continue;
            }

            String proteinId = m.group(1);       // NP_000654 (sin versión)
            char wt = aa3to1(m.group(3));
            int pos = Integer.parseInt(m.group(4));
            char mut = aa3to1(m.group(5));
            if (wt == 0 || mut == 0) {
                writeDiscarded(discardedOut, c);
                continue;
            }

            // Deduplicar por (VariationID, proteinId, wt, pos, mut)
            String key = variationId + "|" + proteinId + "|" + wt + pos + mut;
            if (!seen.add(key)) continue;

            byVariationId.computeIfAbsent(variationId, k -> new ArrayList<>())
                         .add(new ProteinRef(proteinId, wt, pos, mut));
        }

        System.err.println("[HgvsIndex] indexed " + byVariationId.size()
                + " VariationIDs with protein refs");
        return new HgvsIndex(byVariationId);
    }

    private static void writeDiscarded(Writer out, String[] c) throws IOException {
        if (out == null) return;
        out.write(c[COL_SYMBOL]);        out.write('\t');
        out.write(c[COL_VARIATION_ID]);  out.write('\t');
        out.write(c[COL_NUCLEOTIDE]);    out.write('\t');
        out.write(c[COL_PROTEIN_EXPR]);  out.write('\n');
    }

    public List<ProteinRef> findByVariationId(String variationId) {
        if (variationId == null) return Collections.emptyList();
        List<ProteinRef> list = byVariationId.get(variationId);
        return list == null ? Collections.emptyList() : list;
    }

    private static char aa3to1(String aa3) {
        switch (aa3) {
            case "Ala": return 'A'; case "Arg": return 'R'; case "Asn": return 'N';
            case "Asp": return 'D'; case "Cys": return 'C'; case "Gln": return 'Q';
            case "Glu": return 'E'; case "Gly": return 'G'; case "His": return 'H';
            case "Ile": return 'I'; case "Leu": return 'L'; case "Lys": return 'K';
            case "Met": return 'M'; case "Phe": return 'F'; case "Pro": return 'P';
            case "Ser": return 'S'; case "Thr": return 'T'; case "Trp": return 'W';
            case "Tyr": return 'Y'; case "Val": return 'V';
            default: return 0;
        }
    }
}