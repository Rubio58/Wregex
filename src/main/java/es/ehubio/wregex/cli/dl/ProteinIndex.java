package es.ehubio.wregex.cli.dl;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public final class ProteinIndex {

    private static final Pattern GENE = Pattern.compile("GN=(\\S+)");

    private final Map<String, List<ProteinEntry>> byGene;

    private ProteinIndex(Map<String, List<ProteinEntry>> byGene) {
        this.byGene = byGene;
    }

    public static ProteinIndex load(Reader reader) throws IOException {
        Map<String, List<ProteinEntry>> byGene = new HashMap<>();
        BufferedReader br = new BufferedReader(reader);
        StringBuilder seq = new StringBuilder();
        String header = null;
        int count = 0;

        String line;
        while ((line = br.readLine()) != null) {
            if (line.isEmpty()) continue;
            if (line.charAt(0) == '>') {
                if (header != null) { index(header, seq.toString(), byGene); count++; }
                header = line.substring(1);
                seq.setLength(0);
            } else {
                seq.append(line.trim());
            }
        }
        if (header != null) { index(header, seq.toString(), byGene); count++; }

        Comparator<ProteinEntry> order = Comparator
                .comparing((ProteinEntry e) -> !e.isSwissProt())
                .thenComparing(Comparator.comparingInt(ProteinEntry::length).reversed());
        for (List<ProteinEntry> list : byGene.values()) list.sort(order);

        System.err.println("[ProteinIndex] indexed " + count + " proteins, "
                + byGene.size() + " genes");
        return new ProteinIndex(byGene);
    }

    private static void index(String header, String sequence,
                              Map<String, List<ProteinEntry>> byGene) {
        Matcher m = GENE.matcher(header);
        if (!m.find()) return;
        String gene = m.group(1).trim().toUpperCase(Locale.ROOT);

        String[] parts = header.split("\\|", -1);
        String accession = (parts.length >= 3 && !parts[1].isEmpty()) ? parts[1] : null;

        byGene.computeIfAbsent(gene, k -> new ArrayList<>())
              .add(new ProteinEntry(header, sequence, accession));
    }

    public List<ProteinEntry> findByGene(String gene) {
        if (gene == null || gene.isBlank()) return Collections.emptyList();
        List<ProteinEntry> list = byGene.get(gene.trim().toUpperCase(Locale.ROOT));
        return list == null ? Collections.emptyList() : list;
    }
}