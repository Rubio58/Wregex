package es.ehubio.wregex.cli.bl;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import es.ehubio.wregex.cli.dl.ProteinRef;

/**
 * Colapsa VariantResults que difieren solo en proteinId.
 * <p>
 * Criterio estricto: dos filas colapsan si y solo si coinciden en todo
 * excepto el proteinId (y por tanto el change en su parte NP_):
 * <ul>
 *   <li>variant (name, variationId, gene, classification)</li>
 *   <li>status, reason</li>
 *   <li>wtAA + position + mutAA (parte biológica del change)</li>
 *   <li>uniProt, proteinSize</li>
 *   <li>wtWindow, mutWindow</li>
 *   <li>wtMatch, mutMatch, wtScore, mutScore</li>
 *   <li>effect</li>
 * </ul>
 * Las filas colapsadas comparten todos esos campos; solo se diferencian en
 * el proteinId, que se muestra como lista separada por comas.
 */
public final class VariantResultAggregator {

    private VariantResultAggregator() {}

    public static List<VariantResult> collapse(List<VariantResult> input) {
        Map<String, List<VariantResult>> groups = new LinkedHashMap<>();
        for (VariantResult r : input) {
            String key = buildKey(r);
            groups.computeIfAbsent(key, k -> new ArrayList<>()).add(r);
        }

        List<VariantResult> out = new ArrayList<>(groups.size());
        for (List<VariantResult> group : groups.values()) {
            if (group.size() == 1) {
                out.add(group.get(0));
            } else {
                out.add(merge(group));
            }
        }
        return out;
    }

    private static String buildKey(VariantResult r) {
        ProteinRef ref = r.getRef();
        StringBuilder sb = new StringBuilder(256);

        sb.append(safe(r.getVariant().getName())).append('|');
        sb.append(safe(r.getVariant().getVariationId())).append('|');
        sb.append(safe(r.getVariant().getGene())).append('|');
        sb.append(safe(r.getVariant().getClassification())).append('|');
        sb.append(safe(r.getStatus())).append('|');
        sb.append(safe(r.getReason())).append('|');

        if (ref != null) {
            sb.append(ref.getWtAA()).append(ref.getPosition()).append(ref.getMutAA());
        }
        sb.append('|');

        sb.append(safe(r.getUniProt())).append('|');
        sb.append(r.getProteinSize()).append('|');

        sb.append(safe(r.getWtWindow())).append('|');
        sb.append(safe(r.getMutWindow())).append('|');

        sb.append(r.getWtMatch()  == null ? "" : safe(r.getWtMatch().getMatch())).append('|');
        sb.append(r.getMutMatch() == null ? "" : safe(r.getMutMatch().getMatch())).append('|');
        sb.append(r.getWtMatch()  == null ? "" : Double.toString(r.getWtMatch().getScore())).append('|');
        sb.append(r.getMutMatch() == null ? "" : Double.toString(r.getMutMatch().getScore())).append('|');

        sb.append(safe(r.getEffect()));
        return sb.toString();
    }

    private static String safe(String s) {
        return s == null ? "" : s;
    }

    private static VariantResult merge(List<VariantResult> group) {
        VariantResult first = group.get(0);

        List<ProteinRef> refs = new ArrayList<>();
        for (VariantResult r : group)
            refs.addAll(r.getRefs());

        return new VariantResult(
                first.getVariant(),
                refs,
                first.getStatus(),
                first.getReason(),
                first.getWtWindow(),
                first.getMutWindow(),
                first.getWtMatch(),
                first.getMutMatch());
    }
}