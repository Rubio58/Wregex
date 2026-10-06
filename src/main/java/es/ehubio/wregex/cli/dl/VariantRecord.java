package es.ehubio.wregex.cli.dl;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public final class VariantRecord {

    private static final Pattern CHANGE = Pattern.compile("^([A-Z])(\\d+)([A-Z])$");

    private final String name;
    private final List<String>        geneCandidates;
    private final List<ProteinChange> proteinChanges;

    public VariantRecord(String name, String genes, String proteinChange) {
        this.name = name;
        this.geneCandidates = Collections.unmodifiableList(parseMulti(genes));
        this.proteinChanges = Collections.unmodifiableList(parseChanges(proteinChange));
    }

    public String getName() { return name; }
    public List<String> getGeneCandidates() { return geneCandidates; }
    public List<ProteinChange> getProteinChanges() { return proteinChanges; }

    private static List<String> parseMulti(String raw) {
        List<String> out = new ArrayList<>();
        if (raw == null || raw.isBlank()) return out;
        for (String t : raw.split("\\|")) {
            String s = t.trim();
            if (!s.isEmpty()) out.add(s);
        }
        return out;
    }

    private static List<ProteinChange> parseChanges(String raw) {
        List<ProteinChange> out = new ArrayList<>();
        if (raw == null || raw.isBlank()) return out;
        for (String t : raw.split("\\s*,\\s*")) {
            Matcher m = CHANGE.matcher(t.trim());
            if (m.matches())
                out.add(new ProteinChange(
                        m.group(1).charAt(0),
                        Integer.parseInt(m.group(2)),
                        m.group(3).charAt(0)));
        }
        return out;
    }
}