package es.ehubio.wregex.cli.pl;

import java.io.IOException;
import java.io.Writer;
import java.util.List;

import es.ehubio.wregex.cli.bl.VariantResult;

public final class VariantResultWriter {

    private static final String HEADER =
            "variant\tvariationId\tgene\tproteinId\tuniprot\tproteinSize\tchange"
            + "\tstatus\treason\twtWindow\tmutWindow\twtMatch\tmutMatch"
            + "\twtScore\tmutScore\teffect";

    private static final String NA = "N/A";

    private VariantResultWriter() {}

    public static void write(Writer writer, List<VariantResult> results) throws IOException {
        writer.write(HEADER);
        writer.write("\n");
        for (VariantResult r : results) {
            writer.write(toLine(r));
            writer.write("\n");
        }
        writer.flush();
    }

    private static String toLine(VariantResult r) {
        String gene        = r.getVariant().getGene();
        String proteinId   = r.getRef() == null ? null : r.getRef().getProteinId();
        String uniProt     = r.getUniProt();
        String proteinSize = r.getProteinSize() < 0 ? null : String.valueOf(r.getProteinSize());
        String change      = r.getRef() == null ? null : r.getRef().toString();

        return String.join("\t",
                safe(r.getVariant().getName()),
                safe(r.getVariant().getVariationId()),
                safe(gene),
                safe(proteinId),
                safe(uniProt),
                safe(proteinSize),
                safe(change),
                safe(r.getStatus()),
                safe(r.getReason()),
                safe(r.getWtWindow()),
                safe(r.getMutWindow()),
                r.getWtMatch()  == null ? NA : safe(r.getWtMatch().getMatch()),
                r.getMutMatch() == null ? NA : safe(r.getMutMatch().getMatch()),
                r.getWtMatch()  == null ? NA : Double.toString(r.getWtMatch().getScore()),
                r.getMutMatch() == null ? NA : Double.toString(r.getMutMatch().getScore()),
                safe(r.getEffect()));
    }

    private static String safe(String s) {
        return (s == null || s.isEmpty()) ? NA : s;
    }
}