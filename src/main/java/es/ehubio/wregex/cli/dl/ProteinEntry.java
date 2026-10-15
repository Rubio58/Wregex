package es.ehubio.wregex.cli.dl;

/**
 * Entrada de proteína RefSeq. El proteinId se guarda SIN versión.
 * Inmutable.
 */
public final class ProteinEntry {

    private final String proteinId;   // "NP_000654" (sin versión)
    private final String sequence;

    public ProteinEntry(String proteinId, String sequence) {
        this.proteinId = proteinId;
        this.sequence = sequence;
    }

    public String getProteinId() { return proteinId; }
    public String getSequence()  { return sequence; }
    public int    length()       { return sequence.length(); }

    /** Aminoácido en posición 1-based, o '\0' si está fuera de rango. */
    public char aaAt(int position1Based) {
        int i = position1Based - 1;
        if (i < 0 || i >= sequence.length()) return '\0';
        return sequence.charAt(i);
    }
}