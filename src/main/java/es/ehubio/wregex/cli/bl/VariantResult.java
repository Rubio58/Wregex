package es.ehubio.wregex.cli.bl;

import es.ehubio.wregex.Result;
import es.ehubio.wregex.cli.dl.ProteinRef;
import es.ehubio.wregex.cli.dl.VariantRecord;

public final class VariantResult {

    private final VariantRecord variant;
    private final ProteinRef    ref;
    private final String        status;
    private final String        reason;
    private final String        wtWindow;
    private final String        mutWindow;
    private final Result        wtMatch;
    private final Result        mutMatch;

    public VariantResult(VariantRecord variant, ProteinRef ref,
                         String status, String reason,
                         String wtWindow, String mutWindow,
                         Result wtMatch, Result mutMatch) {
        this.variant = variant;
        this.ref = ref;
        this.status = status;
        this.reason = reason;
        this.wtWindow = wtWindow;
        this.mutWindow = mutWindow;
        this.wtMatch = wtMatch;
        this.mutMatch = mutMatch;
    }

    public VariantRecord getVariant() { return variant; }
    public ProteinRef    getRef()     { return ref; }
    public String getStatus()         { return status; }
    public String getReason()         { return reason; }
    public String getWtWindow()       { return wtWindow; }
    public String getMutWindow()      { return mutWindow; }
    public Result getWtMatch()        { return wtMatch; }
    public Result getMutMatch()       { return mutMatch; }

    /** Delega al ref (o null si ref es null). */
    public String getUniProt()     { return ref == null ? null : ref.getUniProt(); }
    public int    getProteinSize() { return ref == null ? -1   : ref.getProteinSize(); }

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