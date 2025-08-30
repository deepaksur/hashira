// Save as ShamirReconstruction.java
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

public class ShamirReconstruction {

    static class Share {
        BigInteger x, y;
        Share(BigInteger x, BigInteger y) { this.x = x; this.y = y; }
    }

    // exact rational using BigInteger
    static class Rational {
        BigInteger num, den;
        Rational(BigInteger n, BigInteger d) {
            if (d.signum() == 0) throw new ArithmeticException("Denominator zero");
            if (d.signum() < 0) { n = n.negate(); d = d.negate(); }
            BigInteger g = n.gcd(d);
            num = n.divide(g); den = d.divide(g);
        }
        static Rational of(BigInteger a) { return new Rational(a, BigInteger.ONE); }
        Rational add(Rational r) { return new Rational(num.multiply(r.den).add(r.num.multiply(den)), den.multiply(r.den)); }
        Rational sub(Rational r) { return new Rational(num.multiply(r.den).subtract(r.num.multiply(den)), den.multiply(r.den)); }
        Rational mul(Rational r) { return new Rational(num.multiply(r.num), den.multiply(r.den)); }
        Rational div(Rational r) { if (r.num.signum()==0) throw new ArithmeticException("Div by zero"); return new Rational(num.multiply(r.den), den.multiply(r.num)); }
        public String toString() { if (den.equals(BigInteger.ONE)) return num.toString(); return num.toString()+"/"+den.toString(); }
        boolean equalsBigInteger(BigInteger b) { return num.equals(b.multiply(den)); }
    }

    // Lagrange interpolation: evaluate polynomial defined by sample at xEval (exact Rational)
    static Rational lagrangeEvalAt(BigInteger xEval, List<Share> sample) {
        Rational total = new Rational(BigInteger.ZERO, BigInteger.ONE);
        int k = sample.size();
        for (int j = 0; j < k; ++j) {
            BigInteger xj = sample.get(j).x;
            BigInteger yj = sample.get(j).y;
            Rational numerator = new Rational(BigInteger.ONE, BigInteger.ONE);
            Rational denom = new Rational(BigInteger.ONE, BigInteger.ONE);
            for (int m = 0; m < k; ++m) {
                if (m == j) continue;
                BigInteger xm = sample.get(m).x;
                numerator = numerator.mul(new Rational(xEval.subtract(xm), BigInteger.ONE));
                denom = denom.mul(new Rational(xj.subtract(xm), BigInteger.ONE));
            }
            Rational ljAt = numerator.div(denom);
            total = total.add(ljAt.mul(Rational.of(yj)));
        }
        return total;
    }

    // Try a subset (indices) as interpolation sample; returns inlier indices or null on failure
    static List<Integer> testSubset(List<Share> allShares, int[] indices) {
        List<Share> sample = new ArrayList<>();
        // ensure unique x in sample
        Set<BigInteger> seen = new HashSet<>();
        for (int idx : indices) {
            BigInteger x = allShares.get(idx).x;
            if (seen.contains(x)) return null;
            seen.add(x);
            sample.add(allShares.get(idx));
        }
        // evaluate polynomial at all x's
        List<Integer> inliers = new ArrayList<>();
        try {
            for (int i = 0; i < allShares.size(); ++i) {
                Share s = allShares.get(i);
                Rational val = lagrangeEvalAt(s.x, sample);
                if (val.equalsBigInteger(s.y)) inliers.add(i);
            }
            return inliers;
        } catch (ArithmeticException e) {
            return null;
        }
    }

    // compute nCk but stop early if greater than threshold
    static long combExceedsThreshold(int n, int k, long threshold) {
        if (k > n) return 0;
        long res = 1;
        for (int i = 1; i <= k; ++i) {
            res = res * (n - k + i) / i;
            if (res > threshold) return res;
        }
        return res;
    }

    // generate combinations (lexicographic) and call consumer for each combination
    static boolean tryAllCombinations(int n, int k, long maxCombs, List<Share> shares, Holder bestHolder) {
        int[] comb = new int[k];
        for (int i = 0; i < k; ++i) comb[i] = i;
        boolean finished = false;
        long tried = 0;
        while (!finished) {
            tried++;
            List<Integer> inliers = testSubset(shares, comb);
            if (inliers != null && inliers.size() > bestHolder.bestInliers.size()) {
                bestHolder.bestInliers = inliers;
                bestHolder.bestSampleIndices = Arrays.copyOf(comb, k);
            }
            // increment combination
            int i = k - 1;
            while (i >= 0 && comb[i] == n - k + i) i--;
            if (i < 0) finished = true;
            else {
                comb[i]++;
                for (int j = i+1; j < k; ++j) comb[j] = comb[j-1] + 1;
            }
            if (tried > maxCombs) break;
        }
        bestHolder.trials = tried;
        return true;
    }

    // random sampling of k-subsets for trials times
    static void randomSampling(int n, int k, int trials, List<Share> shares, Holder bestHolder) {
        Random rnd = new Random(12345);
        for (int t = 0; t < trials; ++t) {
            // sample k distinct indices
            Set<Integer> s = new HashSet<>();
            while (s.size() < k) s.add(rnd.nextInt(n));
            int[] comb = s.stream().mapToInt(Integer::intValue).toArray();
            List<Integer> inliers = testSubset(shares, comb);
            if (inliers != null && inliers.size() > bestHolder.bestInliers.size()) {
                bestHolder.bestInliers = inliers;
                bestHolder.bestSampleIndices = Arrays.copyOf(comb, comb.length);
            }
            bestHolder.trials++;
        }
    }

    static class Holder {
        List<Integer> bestInliers = new ArrayList<>();
        int[] bestSampleIndices = new int[0];
        long trials = 0;
    }

    @SuppressWarnings("unchecked")
    public static void main(String[] args) throws Exception {
        if (args.length < 1) {
            System.out.println("Usage: java -cp .:lib/json-simple-1.1.1.jar ShamirReconstruction input.json");
            return;
        }
        String path = args[0];
        String text = new String(java.nio.file.Files.readAllBytes(java.nio.file.Paths.get(path)), "UTF-8");
        JSONParser p = new JSONParser();
        JSONObject root = (JSONObject)p.parse(text);

        JSONObject keys = (JSONObject)root.get("keys");
        int n = Integer.parseInt(keys.get("n").toString());
        int k = Integer.parseInt(keys.get("k").toString());

        List<Share> shares = new ArrayList<>();
        for (Object objKey : root.keySet()) {
            String sKey = objKey.toString();
            if (sKey.equals("keys")) continue;
            JSONObject entry = (JSONObject)root.get(sKey);
            int base = Integer.parseInt(entry.get("base").toString());
            String value = entry.get("value").toString().trim();
            BigInteger x = new BigInteger(sKey); // key is x coordinate
            BigInteger y = new BigInteger(value, base);
            shares.add(new Share(x, y));
        }

        if (shares.size() != n) {
            System.err.println("Warning: actual parsed shares != n in keys");
        }

        // Decide deterministic vs randomized
        long maxEnumerate = 2000;
        Holder best = new Holder();
        long combCount = combExceedsThreshold(shares.size(), k, maxEnumerate);
        if (combCount <= maxEnumerate) {
            tryAllCombinations(shares.size(), k, maxEnumerate, shares, best);
        } else {
            // randomized
            int trials = 3000;
            randomSampling(shares.size(), k, trials, shares, best);
        }

        if (best.bestSampleIndices.length == 0) {
            System.err.println("Failed to find a good interpolation subset.");
            return;
        }

        // Recompute polynomial using best sample (so we can compute secret and evaluations)
        List<Share> sample = new ArrayList<>();
        for (int idx : best.bestSampleIndices) sample.add(shares.get(idx));
        Rational secret = lagrangeEvalAt(BigInteger.ZERO, sample);

        // evaluate for all shares
        List<Rational> evaluations = new ArrayList<>();
        List<Integer> inliers = new ArrayList<>();
        for (int i = 0; i < shares.size(); ++i) {
            Rational v = lagrangeEvalAt(shares.get(i).x, sample);
            evaluations.add(v);
            if (v.equalsBigInteger(shares.get(i).y)) inliers.add(i);
        }
        List<Integer> outliers = new ArrayList<>();
        for (int i = 0; i < shares.size(); ++i) if (!inliers.contains(i)) outliers.add(i);

        System.out.println("Reconstructed secret f(0) = " + secret.toString());
        System.out.println("Sample used indices (0-based): " + Arrays.toString(best.bestSampleIndices));
        System.out.println("Inliers count: " + inliers.size() + " / " + shares.size());
        System.out.println("Inliers indices: " + inliers);
        System.out.println("Outliers (suspected bad shares): " + outliers);
        System.out.println("Trials performed: " + best.trials);
        System.out.println("\nPer-share evaluation:");
        for (int i = 0; i < shares.size(); ++i) {
            System.out.printf(" %2d: x=%s, y_given=%s, f(x)=%s -> %s\n",
                i, shares.get(i).x.toString(), shares.get(i).y.toString(), evaluations.get(i).toString(),
                inliers.contains(i) ? "OK" : "BAD");
        }
    }
}
