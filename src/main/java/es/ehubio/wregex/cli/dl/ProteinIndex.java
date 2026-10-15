package es.ehubio.wregex.cli.dl;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Reader;
import java.nio.charset.StandardCharsets;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.zip.GZIPInputStream;

/**
 * Índice de proteínas RefSeq, indexadas por proteinID SIN versión.
 * Si varias versiones están presentes, se guarda la de número más alto.
 * <p>
 * Las versiones no se exponen: {@link ProteinEntry} guarda solo el ID sin versión.
 */
public final class ProteinIndex {

    private final Map<String, ProteinEntry> byIdNoVersion;

    private ProteinIndex(Map<String, ProteinEntry> byIdNoVersion) {
        this.byIdNoVersion = byIdNoVersion;
    }

    /**
     * Carga todos los ficheros {@code *.faa.gz} de la carpeta dada, guardando
     * solo las proteínas cuyo ID sin versión esté en {@code wantedIdsNoVersion}.
     */
    public static ProteinIndex loadDirectory(File dir, Set<String> wantedIdsNoVersion) throws IOException {
        File[] files = dir.listFiles((d, n) -> n.endsWith(".faa.gz"));
        if (files == null || files.length == 0)
            throw new IOException("No .faa.gz files in " + dir);

        Map<String, ProteinEntry> map = new HashMap<>();
        Map<String, Integer> bestVersion = new HashMap<>();
        int count = 0;

        for (File f : files) {
            try (Reader r = new InputStreamReader(
                    new GZIPInputStream(new FileInputStream(f)), StandardCharsets.UTF_8)) {
                count += loadInto(r, wantedIdsNoVersion, map, bestVersion);
            }
        }

        System.err.println("[ProteinIndex] indexed " + count + " proteins, "
                + map.size() + " distinct ids");
        return new ProteinIndex(map);
    }

    private static int loadInto(Reader reader, Set<String> wanted,
                                Map<String, ProteinEntry> map,
                                Map<String, Integer> bestVersion) throws IOException {
        BufferedReader br = new BufferedReader(reader);
        StringBuilder seq = new StringBuilder();
        String header = null;
        int count = 0;

        String line;
        while ((line = br.readLine()) != null) {
            if (line.isEmpty()) continue;
            if (line.charAt(0) == '>') {
                if (header != null) {
                    if (store(header, seq.toString(), wanted, map, bestVersion)) count++;
                }
                header = line.substring(1);
                seq.setLength(0);
            } else {
                seq.append(line.trim());
            }
        }
        if (header != null) {
            if (store(header, seq.toString(), wanted, map, bestVersion)) count++;
        }
        return count;
    }

    private static boolean store(String header, String sequence, Set<String> wanted,
                                 Map<String, ProteinEntry> map,
                                 Map<String, Integer> bestVersion) {
        int space = header.indexOf(' ');
        String idWithVersion = space < 0 ? header : header.substring(0, space);
        int dot = idWithVersion.lastIndexOf('.');
        if (dot < 0) return false;

        String idNoVersion = idWithVersion.substring(0, dot);
        if (!wanted.contains(idNoVersion)) return false;

        int version;
        try {
            version = Integer.parseInt(idWithVersion.substring(dot + 1));
        } catch (NumberFormatException e) {
            return false;
        }

        Integer existing = bestVersion.get(idNoVersion);
        if (existing != null && existing >= version) return false;

        bestVersion.put(idNoVersion, version);
        map.put(idNoVersion, new ProteinEntry(idNoVersion, sequence));
        return true;
    }

    /** Devuelve la proteína o null si no está indexada. */
    public ProteinEntry findByProteinId(String proteinId) {
        if (proteinId == null) return null;
        int dot = proteinId.lastIndexOf('.');
        String idNoVersion = dot < 0 ? proteinId : proteinId.substring(0, dot);
        return byIdNoVersion.get(idNoVersion);
    }
}