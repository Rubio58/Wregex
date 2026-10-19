package es.ehubio.wregex.cli.dl;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Reader;
import java.util.ArrayList;
import java.util.List;

public final class VariantTsvReader {

    private static final String SNV        = "single nucleotide variant";
    private static final String MISSENSE   = "missense variant";
    private static final int    MIN_COLUMNS = 16;

    private static final int COL_NAME         = 0;
    private static final int COL_GENES        = 1;
    private static final int COL_VARIATION_ID = 9;
    private static final int COL_TYPE         = 13;
    private static final int COL_CONS         = 14;
    private static final int COL_CLASS        = 15;

    private VariantTsvReader() {}

    public static List<VariantRecord> load(Reader reader) throws IOException {
        List<VariantRecord> records = new ArrayList<>();
        BufferedReader br = new BufferedReader(reader);
        String line;
        boolean headerSkipped = false;

        while ((line = br.readLine()) != null) {
            if (!headerSkipped) { headerSkipped = true; continue; }
            if (line.isBlank()) continue;

            String[] c = line.split("\t", -1);
            if (c.length < MIN_COLUMNS) continue;
            if (!SNV.equalsIgnoreCase(c[COL_TYPE])) continue;
            if (c[COL_CONS] == null || !c[COL_CONS].toLowerCase().contains(MISSENSE)) continue;

            records.add(new VariantRecord(
                    c[COL_NAME],
                    c[COL_GENES],
                    c[COL_VARIATION_ID],
                    c[COL_CLASS]));
        }

        System.err.println("[VariantTsvReader] loaded " + records.size() + " variants");
        return records;
    }
}