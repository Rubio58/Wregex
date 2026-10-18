package es.ehubio.wregex.cli.dl;

/**
 * Una fila del TSV de ClinVar, con los campos necesarios:
 * nombre, VariationID, gen principal y clasificación germinal.
 */
public final class VariantRecord {

    private final String name;
    private final String variationId;
    private final String gene;
    private final String classification;

    public VariantRecord(String name, String genes, String variationId, String classification) {
        this.name = name;
        this.variationId = variationId;
        this.gene = parseMainGene(genes);
        this.classification = classification;
    }

    public String getName() { return name; }
    public String getVariationId() { return variationId; }
    public String getGene() { return gene; }
    public String getClassification() { return classification; }

    private static String parseMainGene(String raw) {
        if (raw == null || raw.isBlank()) return null;
        for (String t : raw.split("\\|")) {
            String s = t.trim();
            if (!s.isEmpty()) return s;
        }
        return null;
    }
}