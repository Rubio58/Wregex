package es.ehubio.wregex.cli.bl;

import es.ehubio.wregex.Result;
import es.ehubio.wregex.cli.dl.ProteinChange;
import es.ehubio.wregex.cli.dl.VariantRecord;

public final class VariantResult {

    private final VariantRecord variant;
    private final ProteinChange change;
    private final String status;    // "RESOLVED" o "SKIPPED"
    private final String reason;    // "" o "NO_PROTEIN" o "WINDOW_INVALID"
    private final String gene;      // gen que resolvió (o vacío)
    private final String uniprot;   // accession (o vacío)
    private final String wtWindow;
    private final String mutWindow;
    private final Result wtMatch;   // null si no hay
    private final Result mutMatch;  // null si no hay

    public VariantResult(VariantRecord variant, ProteinChange change,
                         String status, String reason,
                         String gene, String uniprot,
                         String wtWindow, String mutWindow,
                         Result wtMatch, Result mutMatch) {
        this.variant = variant;
        this.change = change;
        this.status = status;
        this.reason = reason;
        this.gene = gene;
        this.uniprot = uniprot;
        this.wtWindow = wtWindow;
        this.mutWindow = mutWindow;
        this.wtMatch = wtMatch;
        this.mutMatch = mutMatch;
    }

    public VariantRecord getVariant() { return variant; }
    public ProteinChange getChange()  { return change; }
    public String getStatus()         { return status; }
    public String getReason()         { return reason; }
    public String getGene()           { return gene; }
    public String getUniprot()        { return uniprot; }
    public String getWtWindow()       { return wtWindow; }
    public String getMutWindow()      { return mutWindow; }
    public Result getWtMatch()        { return wtMatch; }
    public Result getMutMatch()       { return mutMatch; }

    public String getEffect() {
        if (!"RESOLVED".equals(status)) return "NONE";
        boolean wt = wtMatch != null;
        boolean mut = mutMatch != null;
        if (!wt && mut) return "CREATES";
        if (wt && !mut) return "DESTROYS";
        if (wt && mut) {
            boolean same = wtMatch.getMatch().equals(mutMatch.getMatch())
                    && Double.compare(wtMatch.getScore(), mutMatch.getScore()) == 0;
            return same ? "NONE" : "CHANGES";
        }
        return "NONE";
    }
}