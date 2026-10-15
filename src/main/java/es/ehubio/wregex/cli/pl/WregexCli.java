package es.ehubio.wregex.cli.pl;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Writer;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.List;

import es.ehubio.cli.ArgException;
import es.ehubio.cli.ArgParser;
import es.ehubio.cli.Argument;
import es.ehubio.wregex.Result;
import es.ehubio.wregex.Wregex.WregexException;
import es.ehubio.wregex.cli.bl.SearchOptions;
import es.ehubio.wregex.cli.bl.SearchRunner;
import es.ehubio.wregex.cli.bl.VariantResult;
import es.ehubio.wregex.cli.bl.VariantSearchRunner;
import es.ehubio.wregex.cli.dl.MotifCatalog;
import es.ehubio.wregex.cli.dl.ResultWriter;

public class WregexCli {

    private enum Opt {
        HELP, VERSION, REGEX, PSSM, MOTIF, LIST_MOTIFS, INPUT, OUTPUT,
        THRESHOLD, FILTER_EQUAL, GROUPING, FORMAT
    }

    // ------------------------------------------------------------------
    // main
    // ------------------------------------------------------------------

    public static void main(String[] args) {
        // Dispatch del subcomando "variants"
        if (args.length > 0 && "variants".equals(args[0])) {
            runVariants(Arrays.copyOfRange(args, 1, args.length));
            return;
        }

        ArgParser parser = buildParser();

        if (args.length == 0)
            showUsageAndExit(parser);

        try {
            parser.parseArgs(args);
            if (parser.hasArgument(Opt.HELP.ordinal()))
                showUsageAndExit(parser);
            if (parser.hasArgument(Opt.VERSION.ordinal())) {
                printVersion();
                return;
            }
            if (parser.hasArgument(Opt.LIST_MOTIFS.ordinal())) {
                MotifCatalogPrinter.print(MotifCatalog.load().getMotifs(), System.out);
                return;
            }
            SearchOptions options = toOptions(parser);
            List<Result> results = SearchRunner.run(options);
            try (Writer writer = ResultWriter.open(options.getOutputPath())) {
                ResultWriter.write(options.getFormat(), writer, results);
            }
        } catch (ArgException e) {
            System.err.println(e.getMessage());
            System.err.println();
            System.err.println(parser.getUsage());
            System.exit(1);
        } catch (WregexException e) {
            System.err.println("Invalid configuration: " + e.getMessage());
            System.exit(1);
        } catch (Exception e) {
            System.err.println("Error: " + e.getMessage());
            System.exit(1);
        }
    }

    // ------------------------------------------------------------------
    // Subcomando "variants"
    // ------------------------------------------------------------------

    private static void runVariants(String[] args) {
        ArgParser parser = new ArgParser("wregex variants",
                "Analyze the effect of ClinVar missense pathogenic mutations on a motif");

        Argument help = new Argument(0, 'h', "help", true);
        help.setDescription("Show this help and exit");
        parser.addOption(help);

        Argument motif = new Argument(1, 'm', "motif", true);
        motif.setParamName("MOTIF");
        motif.setDescription("Motif from the built-in catalog (e.g. NES/CRM1, NLS_MonoCore)");
        parser.addOption(motif);

        try {
            parser.parseArgs(args);

            if (parser.hasArgument(0) || args.length == 0) {
                System.out.println(parser.getUsage());
                System.exit(1);
            }
            if (parser.getValue(1) == null) {
                System.err.println("--motif is required");
                System.err.println();
                System.err.println(parser.getUsage());
                System.exit(1);
            }

            String motifToken = parser.getValue(1);

            // Carpeta de salida
            File outputDir = new File("Variants_Results");
            if (!outputDir.exists() && !outputDir.mkdirs()) {
                System.err.println("Could not create directory: " + outputDir.getAbsolutePath());
                System.exit(1);
            }

            // Ejecutar el pipeline
            List<VariantResult> results = new VariantSearchRunner(motifToken, outputDir).run();

            // TSV completo
            String fullName = new File(outputDir, buildOutputFileName(motifToken, "")).getPath();
            try (Writer writer = new FileWriter(fullName)) {
                VariantResultWriter.write(writer, results);
            }
            System.out.println("Output written to " + fullName);

            // TSV solo con effect != NONE
            List<VariantResult> withEffect = new ArrayList<>();
            for (VariantResult r : results)
                if (!"NONE".equals(r.getEffect()))
                    withEffect.add(r);

            String effectName = new File(outputDir, buildOutputFileName(motifToken, "_effects")).getPath();
            try (Writer writer = new FileWriter(effectName)) {
                VariantResultWriter.write(writer, withEffect);
            }
            System.out.println("Effects written to " + effectName
                    + " (" + withEffect.size() + " entries)");

        } catch (ArgException e) {
            System.err.println(e.getMessage());
            System.err.println();
            System.err.println(parser.getUsage());
            System.exit(1);
        } catch (Exception e) {
            System.err.println("Error: " + e.getMessage());
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * Construye un nombre de fichero único a partir del motivo, un sufijo
     * opcional y un timestamp.
     * Ejemplo: motifToken="NES/CRM1", suffix="_effects" -> "NES_CRM1_effects_20261016_101530.tsv".
     */
    private static String buildOutputFileName(String motifToken, String suffix) {
        String sanitized = motifToken.replaceAll("[^A-Za-z0-9]+", "_");
        sanitized = sanitized.replaceAll("^_+|_+$", "");
        if (sanitized.isEmpty())
            sanitized = "motif";

        String timestamp = new SimpleDateFormat("yyyyMMdd_HHmmss").format(new Date());
        return sanitized + suffix + "_" + timestamp + ".tsv";
    }

    // ------------------------------------------------------------------
    // Modo clásico (sin cambios)
    // ------------------------------------------------------------------

    private static void showUsageAndExit(ArgParser parser) {
        System.out.println(parser.getUsage());
        System.exit(1);
    }

    private static void printVersion() {
        try (BufferedReader reader = new BufferedReader(
                new InputStreamReader(WregexCli.class.getResourceAsStream("/version.txt")))) {
            System.out.println("wregex " + reader.readLine());
        } catch (IOException e) {
            System.out.println("wregex (unknown version)");
        }
    }

    private static ArgParser buildParser() {
        ArgParser parser = new ArgParser("wregex", "Protein motif search combining regex and PSSM");

        Argument help = new Argument(Opt.HELP.ordinal(), 'h', "help", true);
        help.setDescription("Show this help and exit");
        parser.addOption(help);

        Argument version = new Argument(Opt.VERSION.ordinal(), 'V', "version", true);
        version.setDescription("Show the application version and exit");
        parser.addOption(version);

        Argument regex = new Argument(Opt.REGEX.ordinal(), 'r', "regex", true);
        regex.setParamName("REGEX");
        regex.setDescription("Wregex regular expression (alternative to --motif)");
        parser.addOption(regex);

        Argument pssm = new Argument(Opt.PSSM.ordinal(), 'p', "pssm", true);
        pssm.setParamName("FILE");
        pssm.setDescription("PSSM file matching the regex capturing groups (only with --regex)");
        parser.addOption(pssm);

        Argument motif = new Argument(Opt.MOTIF.ordinal(), 'm', "motif", true);
        motif.setParamName("MOTIF");
        motif.setDescription("Motif from the built-in catalog, optionally 'name@definition' "
                + "(defaults to the first definition). See --list-motifs for the available names. "
                + "Alternative to --regex/--pssm.");
        parser.addOption(motif);

        Argument listMotifs = new Argument(Opt.LIST_MOTIFS.ordinal(), 'l', "list-motifs", true);
        listMotifs.setDescription("List the built-in motif catalog with its definitions and exit");
        parser.addOption(listMotifs);

        Argument input = new Argument(Opt.INPUT.ordinal(), 'i', "input", true);
        input.setParamName("FILE");
        input.setDescription("Input fasta file");
        parser.addOption(input);

        Argument output = new Argument(Opt.OUTPUT.ordinal(), 'o', "output", true);
        output.setParamName("FILE");
        output.setDescription("Output file (defaults to stdout)");
        parser.addOption(output);

        Argument threshold = new Argument(Opt.THRESHOLD.ordinal(), 't', "threshold", true);
        threshold.setParamName("NUM");
        threshold.setDescription("Minimum PSSM score to keep a match (0 disables filtering)");
        threshold.setDefaultValue(0);
        parser.addOption(threshold);

        Argument filterEqual = new Argument(Opt.FILTER_EQUAL.ordinal(), 'e', "filter-equal", true);
        filterEqual.setDescription("Discard duplicate matches (case-insensitive)");
        parser.addOption(filterEqual);

        Argument grouping = new Argument(Opt.GROUPING.ordinal(), 'g', "grouping", true);
        grouping.setDescription("Collapse overlapping matches, keeping only the best-scoring one");
        parser.addOption(grouping);

        Argument format = new Argument(Opt.FORMAT.ordinal(), 'f', "format", true);
        format.setChoices("csv", "fasta", "aln");
        format.setDescription("Output format. If omitted, it's guessed from --output's extension "
                + "(.fasta/.fa, .aln, otherwise csv); defaults to csv for stdout.");
        parser.addOption(format);

        return parser;
    }

    private static SearchOptions toOptions(ArgParser parser) throws ArgException {
        String regex = parser.getValue(Opt.REGEX.ordinal());
        String pssmPath = parser.getValue(Opt.PSSM.ordinal());
        String motif = parser.getValue(Opt.MOTIF.ordinal());

        if ((regex != null) == (motif != null))
            throw new ArgException("Specify exactly one of --regex or --motif");
        if (motif != null && pssmPath != null)
            throw new ArgException("--pssm cannot be combined with --motif");
        if (parser.getValue(Opt.INPUT.ordinal()) == null)
            throw new ArgException("--input is required");

        String outputPath = parser.getValue(Opt.OUTPUT.ordinal());

        return new SearchOptions(
                regex,
                pssmPath,
                motif,
                parser.getValue(Opt.INPUT.ordinal()),
                outputPath,
                parser.getDoubleValue(Opt.THRESHOLD.ordinal()),
                parser.hasArgument(Opt.FILTER_EQUAL.ordinal()),
                parser.hasArgument(Opt.GROUPING.ordinal()),
                resolveFormat(parser.getValue(Opt.FORMAT.ordinal()), outputPath));
    }

    private static String resolveFormat(String format, String outputPath) {
        if (format != null)
            return format;
        if (outputPath == null)
            return "csv";
        String lower = outputPath.toLowerCase();
        if (lower.endsWith(".fasta") || lower.endsWith(".fa"))
            return "fasta";
        if (lower.endsWith(".aln"))
            return "aln";
        return "csv";
    }
}