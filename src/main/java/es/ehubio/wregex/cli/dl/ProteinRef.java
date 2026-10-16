package es.ehubio.wregex.cli.dl;

/**
 * Referencia a un cambio proteico concreto en una proteína RefSeq.
 * <p>
 * proteinId se guarda SIN versión (p.ej. "NP_000654"). El uniProt y el
 * proteinSize se rellenan después por {@code HgvsIndex.enrich(...)}.
 * <p>
 * Inmutable excepto uniProt y proteinSize, que son package-private setters.
 */
public final class ProteinRef {

    private final String proteinId;   // "NP_000654" (sin versión)
    private final char   wtAA;
    private final int    position;    // 1-based
    private final char   mutAA;

    // Rellenados por HgvsIndex.enrich(...)
    private String uniProt;           // null si no se encontró
    private int    proteinSize;       // -1 si no se encontró

    public ProteinRef(String proteinId, char wtAA, int position, char mutAA) {
        this.proteinId = proteinId;
        this.wtAA = wtAA;
        this.position = position;
        this.mutAA = mutAA;
        this.proteinSize = -1;
    }

    public String getProteinId()   { return proteinId; }
    public char   getWtAA()        { return wtAA; }
    public int    getPosition()    { return position; }
    public char   getMutAA()       { return mutAA; }
    public String getUniProt()     { return uniProt; }
    public int    getProteinSize() { return proteinSize; }

    void setUniProt(String uniProt)      { this.uniProt = uniProt; }
    void setProteinSize(int proteinSize) { this.proteinSize = proteinSize; }

    @Override
    public String toString() {
        return proteinId + ":p." + wtAA + position + mutAA;
    }
}