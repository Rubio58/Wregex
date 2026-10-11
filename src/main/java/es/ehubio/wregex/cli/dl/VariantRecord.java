package es.ehubio.wregex.cli.dl;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public final class VariantRecord {

    /** Notación de 1 letra simple en "Protein change": "Q211E". */
    private static final Pattern CHANGE_1 =
            Pattern.compile("^([A-Z])(\\d+)([A-Z])$");

    /** Notación de 3 letras en el Name: "(p.Pro214Leu)". */
    private static final Pattern NAME_3 =
            Pattern.compile("\\(p\\.([A-Z][a-z]{2})(\\d+)([A-Z][a-z]{2})\\)");

    private static final Map<String, Character> AA3 = Map.ofEntries(
            Map.entry("Ala", 'A'), Map.entry("Arg", 'R'), Map.entry("Asn", 'N'),
            Map.entry("Asp", 'D'), Map.entry("Cys", 'C'), Map.entry("Gln", 'Q'),
            Map.entry("Glu", 'E'), Map.entry("Gly", 'G'), Map.entry("His", 'H'),
            Map.entry("Ile", 'I'), Map.entry("Leu", 'L'), Map.entry("Lys", 'K'),
            Map.entry("Met", 'M'), Map.entry("Phe", 'F'), Map.entry("Pro", 'P'),
            Map.entry("Ser", 'S'), Map.entry("Thr", 'T'), Map.entry("Trp", 'W'),
            Map.entry("Tyr", 'Y'), Map.entry("Val", 'V'));

    private final String name;
    private final List<String>        geneCandidates;
    private final List<ProteinChange> proteinChanges;

    public VariantRecord(String name, String genes, String proteinChange) {
        this.name = name;
        this.geneCandidates = Collections.unmodifiableList(parseMulti(genes));
        this.proteinChanges = Collections.unmodifiableList(parseChanges(proteinChange, name));
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

    /**
     * Parsea los cambios proteicos. Orden de preferencia:
     * <ol>
     *   <li>"Protein change" (notación 1 letra, posible lista "Q160R, Q193R").</li>
     *   <li>Si lo anterior no produce resultados, "Name" (notación 3 letras
     *       {@code (p.Pro214Leu)}, un solo cambio).</li>
     * </ol>
     * Siempre devuelve {@link ProteinChange} en notación de 1 letra.
     * Si ninguna fuente produce resultados, devuelve lista vacía.
     */
    private static List<ProteinChange> parseChanges(String proteinChangeCol, String name) {
        // 1) "Protein change" primero
        if (proteinChangeCol != null && !proteinChangeCol.isBlank()) {
            List<ProteinChange> out = new ArrayList<>();
            for (String token : proteinChangeCol.split("\\s*,\\s*")) {
                Matcher m = CHANGE_1.matcher(token.trim());
                if (m.matches()) {
                    out.add(new ProteinChange(
                            m.group(1).charAt(0),
                            Integer.parseInt(m.group(2)),
                            m.group(3).charAt(0)));
                }
            }
            if (!out.isEmpty()) return out;
        }

        // 2) Fallback: "Name" (notación 3 letras)
        if (name != null) {
            Matcher m = NAME_3.matcher(name);
            if (m.find()) {
                Character wt = AA3.get(m.group(1));
                Character mut = AA3.get(m.group(3));
                if (wt != null && mut != null) {
                    List<ProteinChange> out = new ArrayList<>(1);
                    out.add(new ProteinChange(wt, Integer.parseInt(m.group(2)), mut));
                    return out;
                }
            }
        }

        // 3) Nada
        return Collections.emptyList();
    }
}