package es.ehubio.wregex.cli.dl;

/**
 * Referencia a un cambio proteico concreto en una proteína RefSeq.
 * <p>
 * El proteinId se guarda SIN versión (p.ej. "NP_000654"), porque la versión
 * solo se resuelve al buscar en el FASTA. La versión más alta encontrada
 * es la que se usa para extraer la ventana.
 * <p>
 * Inmutable.
 */
public final class ProteinRef {

    private final String proteinId;   // "NP_000654" (sin versión)
    private final char   wtAA;
    private final int    position;    // 1-based
    private final char   mutAA;

    public ProteinRef(String proteinId, char wtAA, int position, char mutAA) {
        this.proteinId = proteinId;
        this.wtAA = wtAA;
        this.position = position;
        this.mutAA = mutAA;
    }

    public String getProteinId() { return proteinId; }
    public char   getWtAA()      { return wtAA; }
    public int    getPosition()  { return position; }
    public char   getMutAA()     { return mutAA; }

    @Override
    public String toString() {
        return proteinId + ":p." + wtAA + position + mutAA;
    }
}