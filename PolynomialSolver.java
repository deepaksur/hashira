import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.math.BigInteger;
import java.util.*;
import org.json.JSONObject;

public class PolynomialSolver {

    // Convert value in given base to decimal BigInteger
    static BigInteger convertBase(String value, int base) {
        return new BigInteger(value, base);
    }

    // Lagrange Interpolation to get polynomial coefficients
    static BigInteger[] lagrangeInterpolation(int[] x, BigInteger[] y, int k) {
        BigInteger[] coeffs = new BigInteger[k];
        Arrays.fill(coeffs, BigInteger.ZERO);

        for (int i = 0; i < k; i++) {
            // Build Lagrange basis polynomial L_i(x)
            BigInteger denom = BigInteger.ONE;
            BigInteger[] term = new BigInteger[k];
            Arrays.fill(term, BigInteger.ZERO);
            term[0] = BigInteger.ONE;

            for (int j = 0; j < k; j++) {
                if (i != j) {
                    denom = denom.multiply(BigInteger.valueOf(x[i] - x[j]));

                    // Multiply polynomial by (X - xj)
                    BigInteger[] newTerm = new BigInteger[k];
                    Arrays.fill(newTerm, BigInteger.ZERO);

                    for (int d = k - 2; d >= 0; d--) {
                        if (term[d] != null) {
                            newTerm[d + 1] = newTerm[d + 1].add(term[d]); // coeff for x^(d+1)
                            newTerm[d] = newTerm[d].subtract(term[d].multiply(BigInteger.valueOf(x[j])));
                        }
                    }
                    term = newTerm;
                }
            }

            // Multiply L_i(x) by y[i]/denom
            for (int d = 0; d < k; d++) {
                coeffs[d] = coeffs[d].add(y[i].multiply(term[d]).divide(denom));
            }
        }
        return coeffs;
    }

    public static void main(String[] args) throws IOException {
        // Read JSON from file
        String jsonInput = new String(Files.readAllBytes(Paths.get("input.json")));
        JSONObject obj = new JSONObject(jsonInput);

        int n = obj.getJSONObject("keys").getInt("n");
        int k = obj.getJSONObject("keys").getInt("k");

        int[] x = new int[n];
        BigInteger[] y = new BigInteger[n];
        int idx = 0;

        for (String key : obj.keySet()) {
            if (!key.equals("keys")) {
                int xi = Integer.parseInt(key);
                JSONObject entry = obj.getJSONObject(key);
                int base = entry.getInt("base");
                String value = entry.getString("value");
                BigInteger yi = convertBase(value, base);

                x[idx] = xi;
                y[idx] = yi;
                idx++;
            }
        }

        // Take first k points (or any valid k points)
        int[] xk = Arrays.copyOf(x, k);
        BigInteger[] yk = Arrays.copyOf(y, k);

        BigInteger[] coeffs = lagrangeInterpolation(xk, yk, k);

        // Print coefficients
        System.out.println("Polynomial Coefficients:");
        for (int i = 0; i < coeffs.length; i++) {
            System.out.println("a" + i + " = " + coeffs[i]);
        }
    }
}