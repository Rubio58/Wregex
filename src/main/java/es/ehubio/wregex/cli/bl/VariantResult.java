package es.ehubio.wregex.cli.bl;

import java.util.Collections;
import java.util.List;

import es.ehubio.wregex.Result;
import es.ehubio.wregex.cli.dl.ProteinRef;
import es.ehubio.wregex.cli.dl.VariantRecord;

public final class VariantResult {

    private final VariantRecord    variant;
    private final List<ProteinRef> refs;
    private final String           status;
    private final String           reason;
    private final String           wtWindow;
    private final String           mutWindow;
    private final Result           wtMatch;
    private final Result           mutMatch;

    public VariantResult(VariantRecord variant, List<ProteinRef> refs,
                         String status, String reason,
                         String wtWindow, String mutWindow,
                         Result wtMatch, Result mutMatch) {
        this.variant = variant;
        this.refs = refs == null ? Collections.emptyList() : Collections.unmodifiableList(refs);
        this.status = status;
        this.reason = reason;
        this.wtWindow = wtWindow;
        this.mutWindow = mutWindow;
        this.wtMatch = wtMatch;
        this.mutMatch = mutMatch;
    }

    public VariantRecord    getVariant() { return variant; }
    public List<ProteinRef> getRefs()    { return refs; }
    public ProteinRef       getRef()     { return refs.isEmpty() ? null : refs.get(0); }
    public String getStatus()            { return status; }
    public String getReason()            { return reason; }
    public String getWtWindow()          { return wtWindow; }
    public String getMutWindow()         { return mutWindow; }
    public Result getWtMatch()           { return wtMatch; }
    public Result getMutMatch()          { return mutMatch; }

    public String getUniProt() {
        for (ProteinRef r : refs)
            if (r.getUniProt() != null) return r.getUniProt();
        return null;
    }

    public int getProteinSize() {
        for (ProteinRef r : refs)
            if (r.getProteinSize() >= 0) return r.getProteinSize();
        return -1;
    }

    public String getEffect() {
        if (!"RESOLVED".equals(status)) return "NONE";
        boolean wt = wtMatch != null;
        boolean mut = mutMatch != null;
        if (!wt && mut) return "CREATES";
        if (wt && !mut) return "DESTROYS";
        if (wt && mut) {
            boolean sameScore = Double.compare(wtMatch.getScore(), mutMatch.getScore()) == 0;
            return sameScore ? "NONE" : "CHANGES";
        }
        return "NONE";
    }
}