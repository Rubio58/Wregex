package es.ehubio.wregex.cli.dl;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Reader;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

/**
 * Mapa de proteinIDs RefSeq (NP_XXXXXX, sin versión) a accesiones UniProt.
 * <p>
 * Fuente: HUMAN_9606_idmapping_selected.tab.gz.
 * Columnas (0-based): 0 = UniProtKB-AC, 3 = RefSeq (NP_...; NP_...).
 * Si un NP_ aparece en varias entradas UniProt, se guarda el primero.
 */
public final class UniProtMap {

    private static final int COL_UNIPROT = 0;
    private static final int COL_REFSEQ  = 3;

    private final Map<String, String> byProteinId;

    private UniProtMap(Map<String, String> byProteinId) {
        this.byProteinId = byProteinId;
    }

    public static UniProtMap load(Reader reader, Set<String> wanted) throws IOException {
        Map<String, String> map = new HashMap<>();
        BufferedReader br = new BufferedReader(reader);
        String line;

        while ((line = br.readLine()) != null) {
            if (line.isEmpty()) continue;
            String[] c = line.split("\t", -1);
            if (c.length < 4) continue;

            String uniprot = c[COL_UNIPROT];
            String refseqField = c[COL_REFSEQ];
            if (uniprot.isEmpty() || refseqField.isEmpty() || refseqField.equals("-"))
                continue;

            for (String token : refseqField.split(";")) {
                String np = token.trim();
                if (np.isEmpty()) continue;
                int dot = np.lastIndexOf('.');
                String npNoVersion = dot < 0 ? np : np.substring(0, dot);
                if (!wanted.contains(npNoVersion)) continue;
                map.putIfAbsent(npNoVersion, uniprot);
            }
        }

        System.err.println("[UniProtMap] mapped " + map.size() + " protein ids to uniprot");
        return new UniProtMap(map);
    }

    /** Devuelve el UniProt accession o null si no hay mapeo. */
    public String get(String proteinIdNoVersion) {
        return proteinIdNoVersion == null ? null : byProteinId.get(proteinIdNoVersion);
    }
}