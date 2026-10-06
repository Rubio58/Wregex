package es.ehubio.wregex.cli.dl;

public final class ProteinEntry {

    private final String header;
    private final String sequence;
    private final String accession;

    public ProteinEntry(String header, String sequence, String accession) {
        this.header = header;
        this.sequence = sequence;
        this.accession = accession;
    }

    public String getSequence()  { return sequence; }
    public String getAccession() { return accession; }
    public int    length()       { return sequence.length(); }

    public char aaAt(int position1Based) {
        int i = position1Based - 1;
        if (i < 0 || i >= sequence.length()) return '\0';
        return sequence.charAt(i);
    }

    boolean isSwissProt() { return header.startsWith("sp|"); }
}