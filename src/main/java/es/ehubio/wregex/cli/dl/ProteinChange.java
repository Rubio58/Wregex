package es.ehubio.wregex.cli.dl;

public final class ProteinChange {
    private final char wtAA;
    private final int  position;   // 1-based
    private final char mutAA;

    public ProteinChange(char wtAA, int position, char mutAA) {
        this.wtAA = wtAA;
        this.position = position;
        this.mutAA = mutAA;
    }

    public char getWtAA()     { return wtAA; }
    public int  getPosition() { return position; }
    public char getMutAA()    { return mutAA; }

    @Override
    public String toString() { return "" + wtAA + position + mutAA; }
}