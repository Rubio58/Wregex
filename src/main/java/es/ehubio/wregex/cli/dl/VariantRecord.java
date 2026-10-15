package es.ehubio.wregex.cli.dl;

/**
 * Una fila del TSV de ClinVar, con los campos necesarios:
 * nombre, VariationID y gen principal.
 */
public final class VariantRecord {

    private final String name;
    private final String variationId;
    private final String gene;   // gen principal, o null si no hay

    public VariantRecord(String name, String genes, String variationId) {
        this.name = name;
        this.variationId = variationId;
        this.gene = parseMainGene(genes);
    }

    public String getName() { return name; }
    public String getVariationId() { return variationId; }
    public String getGene() { return gene; }

    /**
     * Extrae el primer gen no vacío de la columna "Gene(s)".
     * Ejemplo: "AAAS|LOC130007973" -> "AAAS".
     */
    private static String parseMainGene(String raw) {
        if (raw == null || raw.isBlank()) return null;
        for (String t : raw.split("\\|")) {
            String s = t.trim();
            if (!s.isEmpty()) return s;
        }
        return null;
    }
}