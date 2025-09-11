import java.io.IOException;
import java.math.BigInteger;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class MyClass {

    // ---------- Fraction ----------
    private static final class Frac {
        final BigInteger n;
        final BigInteger d;
        Frac(BigInteger n, BigInteger d) {
            if (d.signum() == 0) throw new IllegalArgumentException("Zero denominator");
            if (d.signum() < 0) { n = n.negate(); d = d.negate(); }
            BigInteger g = n.gcd(d);
            this.n = n.divide(g);
            this.d = d.divide(g);
        }
        Frac add(Frac o) {
            return new Frac(this.n.multiply(o.d).add(o.n.multiply(this.d)), this.d.multiply(o.d));
        }
        Frac mul(Frac o) { return new Frac(this.n.multiply(o.n), this.d.multiply(o.d)); }
        Frac div(Frac o) { return new Frac(this.n.multiply(o.d), this.d.multiply(o.n)); }
        BigInteger toBigIntExact() {
            if (!d.equals(BigInteger.ONE)) throw new ArithmeticException("Non-integer result");
            return n;
        }
        static Frac of(BigInteger x) { return new Frac(x, BigInteger.ONE); }
    }

    private static final class Point {
        final BigInteger x;
        final BigInteger y;
        Point(BigInteger x, BigInteger y) { this.x = x; this.y = y; }
    }

    private static class TestCase {
        final int n;
        final int k;
        final List<Point> pts;
        TestCase(int n, int k, List<Point> pts) { this.n = n; this.k = k; this.pts = pts; }
    }

    private static TestCase parseTestCase(String json) {
        Pattern keysPat = Pattern.compile(
            "\"keys\"\\s*:\\s*\\{[^}]*?\"n\"\\s*:\\s*(\\d+)\\s*,\\s*\"k\"\\s*:\\s*(\\d+)[^}]*\\}",
            Pattern.DOTALL);
        Matcher km = keysPat.matcher(json);
        if (!km.find()) throw new IllegalArgumentException("Cannot find n/k");
        int n = Integer.parseInt(km.group(1));
        int k = Integer.parseInt(km.group(2));

        Pattern entryPat = Pattern.compile(
            "\"(\\d+)\"\\s*:\\s*\\{\\s*\"base\"\\s*:\\s*\"([0-9A-Za-z]+)\"\\s*,\\s*\"value\"\\s*:\\s*\"([0-9A-Za-z]+)\"\\s*\\}",
            Pattern.DOTALL);
        Matcher em = entryPat.matcher(json);

        List<Point> pts = new ArrayList<>();
        while (em.find()) {
            BigInteger x = new BigInteger(em.group(1));
            int radix = Integer.parseInt(em.group(2));
            BigInteger y = new BigInteger(em.group(3), radix);
            pts.add(new Point(x, y));
        }
        pts.sort(Comparator.comparing(p -> p.x));
        return new TestCase(n, k, pts);
    }

    private static BigInteger constantTermAtZero(List<Point> points, int k) {
        List<Point> subset = points.subList(0, k);
        Frac c = Frac.of(BigInteger.ZERO);
        for (int i = 0; i < subset.size(); i++) {
            BigInteger xi = subset.get(i).x;
            BigInteger yi = subset.get(i).y;
            Frac num = Frac.of(BigInteger.ONE), den = Frac.of(BigInteger.ONE);
            for (int j = 0; j < subset.size(); j++) {
                if (i == j) continue;
                num = num.mul(Frac.of(subset.get(j).x.negate()));
                den = den.mul(Frac.of(xi.subtract(subset.get(j).x)));
            }
            Frac term = Frac.of(yi).mul(num.div(den));
            c = c.add(term);
        }
        return c.toBigIntExact();
    }

    public static void main(String[] args) {
        try {
            // read sample.json from current folder
            String json = Files.readString(Paths.get("sample.json"));
            TestCase tc = parseTestCase(json);
            BigInteger secret = constantTermAtZero(tc.pts, tc.k);
            System.out.println("Secret c = " + secret);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
