using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Security;
using System.Text;

// Written by CorruptConverter, discord: corruptconverter#0000
public class EntryPoint
{
    // This is where the code starts running
    public static void Main()
    {
        try
        {
            Console.OutputEncoding = Encoding.UTF8;
        }
        catch
        {
        }

        // ##### Enter code starting here
        Vect.DefaultSigfigs = 8;

        uint genealogy = (uint)Global.rand.Next(-2147483648, 2147483647);

        Tessarine expTestCmp = Tessarine.Random(8, -10, 10);
        //Console.WriteLine(new Tessarine(0, 0.1234567890123456789, 1.234567890123456789, 12.34567890123456789, 123.4567890123456789, 1234.567890123456789, 12345.67890123456789, 123456.7890123456789, 1234567.890123456789, 12345678.90123456789, 123456789.0123456789, 1234567890.123456789));
        //Console.WriteLine("Compare: "+(TessMath.exp(expTestCmp)-TessMath.expTest(expTestCmp)).AbsMagnitude());
        Console.WriteLine("Exp:" + TessMath.exp(expTestCmp));
        Console.WriteLine();
        Console.WriteLine("ExpTest: " + (TessMath.expTest(expTestCmp)));
        Console.WriteLine("PerfTest: Exp");
        {
            double time = 0;
            for (int i = 0; i < 1; i++) // Set to 1 due to dotnetfiddle limitations. Likely runs much better locally.
            {
                /*for (int j = 0; j < 4096; j++)
                {
                    Math.Cosh(i*j);
                    Math.Sinh(i*j);
                }
                goto Pen;*/
                Tessarine tessTest = Tessarine.Random(4096, -10, 10);
                Stopwatch stp = new Stopwatch();
                stp.Start();
                TessMath.expTest(tessTest);
                stp.Stop();
                time += stp.ElapsedMilliseconds;
            Pen:
                Console.WriteLine("Finished: " + i);
            }
            Console.WriteLine("Average time: " + (time / 1000));
            Console.WriteLine("Total time: " + time);
        }

        return;

        Console.WriteLine("Started");
        Console.WriteLine("Test: " + new HybridTessarine(13, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16).Abs());
        HybridTessarine t = HybridTessarine.Random(genealogy, 16, -10, 10);  // Creates a random 256-dimensional Tessarine with each component ranging from -10 to +10
        HybridTessarine t2 = HybridTessarine.Random(genealogy, 16, -10, 10); // Creates another random 256-dimensional Tessarine
        Tessarine tP = t;
        Tessarine t2P = t2;
        Tessarine prod = tP * t2P;
        Console.WriteLine("Diff: " + ((t * t2) - new HybridTessarine(genealogy, prod)));
        Console.WriteLine("Diff to/from: " + (t - new HybridTessarine(genealogy, tP)));
        Console.WriteLine("Genealogy: " + Convert.ToString(genealogy, 2));
        Console.WriteLine("Conjugate product: " + (t * t.conj()));
        Console.WriteLine("Value 1: " + t);
        Console.WriteLine("Value 2: " + t2);
        Console.WriteLine("t*t2 = " + (t * t2) + "\n");

        Console.WriteLine("Proper[t] = " + tP);
        Console.WriteLine("Proper[t2] = " + t2P);
        Console.WriteLine("tP * t2P = " + prod);

        Console.WriteLine("Hybrid[prod] = " + new HybridTessarine(genealogy, prod));

        Console.WriteLine("Difference in multiplication: " + ((t * t2) - new HybridTessarine(genealogy, (((Tessarine)t) * ((Tessarine)t2)))).AbsMagnitude());
        Console.WriteLine("Absolute values difference: " + (t.Abs() - tP.Abs()));
        Console.WriteLine("Ended");

        //Tessarine lw = TessMath.LambertW(t);
        //Console.WriteLine("Error in LambertW(t) = " + ((lw * TessMath.exp(lw)) - t).AbsMagnitude());

        double absVal = t.Abs();
        double absVal2 = t2.Abs();

        Console.WriteLine();

        Console.WriteLine("|t| = " + absVal);
        Console.WriteLine("|t2| = " + absVal2);
        Console.WriteLine("|t|*|t2| = " + (absVal * absVal2));
        Console.WriteLine("|t*t2| = " + (t * t2).Abs());

        Console.WriteLine("Abs test 1:");

        Console.WriteLine("|3.124121 + 5.124243i -3.656544j + 5.222234k| = " + new Tessarine(3.124121, 5.124243, -3.656544, 5.222234).Abs());

        Console.WriteLine();

        Console.WriteLine("Abs Test 2: ");
        double absValTest = TessMath.AbsTest(t);
        double absVal2Test = TessMath.AbsTest(t2);

        Console.WriteLine("|t| = " + absValTest);
        Console.WriteLine("|t2| = " + absVal2Test);
        Console.WriteLine("|t|*|t2| = " + (absValTest * absVal2Test));
        Console.WriteLine("|t*t2| = " + TessMath.AbsTest(t * t2));

        Tessarine logResult = TessMath.log(t, t2); // Saves the result of log_t(t2) to a variable called "logResult"
        Tessarine powResult = TessMath.pow(t, logResult); // Saves the result of t^logResult to a variable called "powResult"
        Tessarine diff = powResult - t2; // Computes the difference and saves it in "diff"
        double diffAM = diff.AbsMagnitude();
        // Computes the sum of absolute values of all vectors in diff (Giving a good indicator of the error in our calculations)

        Console.WriteLine("t = " + t); // Outputs Tessarine 't' to the console
        Console.WriteLine("t2 = " + t2); // Outputs Tessarine 't2' to the console
        Console.WriteLine("Error in t^(log_t(t2))=" + diffAM); // Outputs the result diffAM to the console below when you press the run button

        /*
            Result might be in scientific notation if the result is either large enough or small enough, an example is:
		    "3.66918474770106E-12"

		    This is just 3.66918474770106 * 10^(-12)
		
		    By default, I set the tessarines to round to 9 significant figures, to get less numbers that look something like 5.99999999999999873
            as a result of floating point error.
        
            To change this, add the following line somewhere here (Not in a comment):
            Vect.DefaultSigfigs = 13;

            And replace '13' with however many significant figures you want.
        */

        // <-- This is a single-line comment; It lets you write notes for the reader, rather than code for the computer

        /*
			This is a multi-line comment
			In general, if it's green, it's probably a comment.
		*/

        // ##### Enter code ending here
    }
}
/*
    Available functions:
    TessMath.exp(a) . . . .>> e^a
    TessMath.ln(a)  . . . .>> Natural log of 'a'
    TessMath.log(b,a) . . .>> Log base 'b' of 'a' ['b' comes before 'a']
    TessMath.atan2(a,b) . .>> 2-argument arctan
    TessMath.conj(a)  . . .>> Conjugate of 'a'
    TessMath.div(a,b) . . .>> a/b
    TessMath.sqrt(a)  . . .>> Square-root of 'a'
    TessMath.pow(a,b) . . .>> a^b

    TessMath.cos(a)
    TessMath.sin(a)
    TessMath.tan(a)
    TessMath.cot(a)
    TessMath.sec(a)
    TessMath.csc(a)

    TessMath.acos(a)
    TessMath.asin(a)
    TessMath.atan(a)
    TessMath.acot(a)
    TessMath.asec(a)
    TessMath.acsc(a)

    TessMath.cosh(a)
    TessMath.sinh(a)
    TessMath.tanh(a)
    TessMath.coth(a)
    TessMath.sech(a)
    TessMath.csch(a)

    TessMath.acosh(a)
    TessMath.asinh(a)
    TessMath.atanh(a)
    TessMath.acoth(a)
    TessMath.asech(a)
    TessMath.acsch(a)

    Define tessarine as:
    Tessarine variableName = new Tessarine(a,b,c,d,e,f, ...);
    // Parameters are all doubles (Or all Vects)

    variableName.AbsMagnitude() -> Sum of absolute values of every vector in 'variableName'

    variableName.imul() -> Multiply variableName by 'i' (First imaginary unit)

    variableName.left -> Left half of Tessarine 'variableName'
    variableName.right -> Right half of Tessarine 'variableName'

    variableName.IsZero() -> Returns true if every vector in 'variableName' is zero

    Tessarines can be compared by a==b and a!=b,
	multiplied by a*b, divided by a/b,
	added together by a+b, subtracted by a-b,
	and exponentiated with a^b.
	Interoperability between Tessarines, doubles, and Vects is supported to some extent.
*/

// Library of functions that you can run on Tessarines and shit
public static class TessMath
{
    public static double AbsTest(Tessarine t)
    {
        if (t.dimensions == 2)
        {
            return Math.Sqrt(t[0].value * t[0].value + t[1].value * t[1].value);
        }
        return Math.Sqrt(AbsTest(t.left + t.right) * AbsTest(t.left - t.right));
    }
    // Compute the exp (e^x) of a single Vector
    public static Tessarine exp(Vect v)
    {
        if (double.IsNaN(v.value))
        {
            return Tessarine.NaN;
        }
        return new Tessarine(expSplit(v));
    }
    // Returns the two parts of the exp (e^x) of a single vector; The real and (usually) imaginary part
    public static Vect[] expSplit(Vect v)
    {
        if (v.value == 0)
        {
            return new Vect[] { new Vect(1, 0), new Vect(0, 1) };
        }
        Vect c, s;
        if ((v.index & 1) == 0)
        {
            c = new Vect(Math.Cosh(v.value), 0);
            s = new Vect(Math.Sinh(v.value), v.index);
        }
        else
        {
            c = new Vect(Math.Cos(v.value), 0);
            s = new Vect(Math.Sin(v.value), v.index);
        }
        return new Vect[] { c, s };
    }
    public static void expSplit(Vect v, out double real, out double imaginary)
    {
        if (v.value == 0)
        {
            real = 1;
            imaginary = 0;
        }
        if (v.index == 0)
        {
            real = Math.Exp(v.value);
            imaginary = new Vect(0, 1);
        }
        if ((v.index & 1) == 0)
        {
            real = Math.Cosh(v.value);
            imaginary = Math.Sinh(v.value);
        }
        else
        {
            real = Math.Cos(v.value);
            imaginary = Math.Sin(v.value);
        }
    }

    public static Tessarine expTest(Tessarine t)
    {
        if (t.IsNaN())
        {
            return Tessarine.NaN;
        }
        Vect[] v = t.vectors;
        double[] res = new double[v.Length],
                 temp1 = new double[v.Length],
                 temp2;
        //res[0] = Math.Exp(v[0].value);
        res[0] = 1;

        for (int i = 0; i < v.Length; i++)
        {
            double vM0;
            double vM1;
            TessMath.expSplit(v[i], out vM0, out vM1);

            int len = i.CeilDimension();

            for (int j = 0; j < len; j++)
            {
                temp1[j] = res[j] * vM0;
            }
            for (int j = 0; j < len; j++)
            {
                if ((j & i & 1) == 0)
                {
                    temp1[j ^ i] += (res[j] * vM1);
                }
                else
                {
                    temp1[j ^ i] -= (res[j] * vM1);
                }
            }
            temp2 = res;
            res = temp1;
            temp1 = temp2;
        }
        return new Tessarine(res, true);
    }

    // Compute the exp (e^x) of a Tessarine
    public static Tessarine exp(Tessarine t)
    {
        if (t.IsNaN())
        {
            return Tessarine.NaN;
        }
        Vect[] v = t.vectors;
        Tessarine res = new Tessarine(new double[t.dimensions]);
        res[0] = new Vect(1, 0);

        for (int i = 0; i < v.Length; i++)
        {
            res *= TessMath.expSplit(v[i]);
        }
        return res;
    }
    public static HybridTessarine exp(HybridTessarine t)
    {
        return new HybridTessarine(t.genealogy, exp((Tessarine)t));
    }

    // Precomputed constants to speed up computation
    private const double halfpi = 1.57079632679489661923132169163975;
    private const double rescaleFittingA = 512;
    private const double rescaleFittingB = 0.001953125;
    private const double rescaleFittingC = 6.2383246250395077847550890931235891126795;

    //private const double cos1 = 0.5403023058681397174009366074429766;
    //private const double sin1 = 0.841470984807896506652502321630298999622563;

    // Compute the natural log of a Tessarine
    public static Tessarine ln(Tessarine t)
    {
        return LogRecursive(t, -1);
    }
    public static HybridTessarine ln(HybridTessarine t)
    {
        return new HybridTessarine(t.genealogy, ln((Tessarine)t));
    }

    // Internal recursive implementation of Natural Log
    private static Tessarine LogRecursive(Tessarine t, int debug = 0)
    {
        if (t.IsNaN())
        {
            return Tessarine.NaN;
        }
        t = t.Reduce(); // Reduce the tessarine to a lower-dimension, if possible
        // This line only begins to do anything if the right half of the tessarine is all zeroes
        if (t.IsZero()) // Handle the special case that the tessarine is zero.
        {
            // If I figure out how to kind-of handle when a tessarine is an impossible, a+aj or a-aj, I'll add it
            return new Tessarine(double.NegativeInfinity);
        }
        if (t.IsZeroDivisor())
        {
            return Tessarine.NaN;
        }
        // Handle complex natural logarithm in a specific way;
        // This is the entire reason why this function returns a result at all in the end
        if (t.dimensions == 2)
        {
            double a1 = t[0].value;
            double b1 = t[1].value;
            return new Tessarine(Math.Log(Math.Sqrt(a1 * a1 + b1 * b1)), Math.Atan2(b1, a1));
            // ln(|t|)+i*arg(t)
        }
        double absMag = t.AbsMagnitude(); // Guesstimate the magnitude of 't' by summing the absolute-values of each term
        double shiftLogResult = 0; // This will be added to the result in the end to reverse the shifting that may happen next
        while (absMag > rescaleFittingA) // If it's too big, make it smaller by t/512
        {
            t = t * rescaleFittingB;
            absMag *= rescaleFittingB;
            shiftLogResult += rescaleFittingC;
        }
        while (absMag < rescaleFittingB) // If it's too small, make it bigger by t*512
        {
            t = t * rescaleFittingA;
            absMag *= rescaleFittingA;
            shiftLogResult -= rescaleFittingC;
        }
        /*
            All the notes I had written to derive a better natural-log that worked on large tessarines when using doubles.
			If someone *really* wants to know, maybe I'll explain it to them. If it's Ed, I'll for sure go through the full explanation.
			Slycedf? 2yce? Jerridium? Idk, maybe I'll explain it to you.

            X=ln(b-a)
            P=ln(b+a)

            ln(q=b+aj)
            # z=-i*atan2(ai,b)
            = ln(b/cosh(z))+jz
            = ln(b)-ln(cosh(z))+jz

            atan2(a,b)
            = -i*ln([b+ia]/sqrt(a^2+b^2))

            -i*atan2(a,b)
            = -ln([b+ia]/sqrt(a^2 + b^2))

            -i*atan2(ia,b)
            = -ln([b-a]/sqrt(b^2 - a^2))
            = -(ln(b-a)-0.5*ln(b^2 - a^2))
            = -(X-0.5(X+P))
            = -(X-0.5X-0.5P)
            = -0.5(X-P) = -i*atan2(ia,b)

            ln[  (b-a)/sqrt(b^2 - a^2)  ]
            
            ln(sqrt(b^2 - a^2))
            = 0.5*ln(b^2 - a^2)
            = 0.5*(ln(b-a) + ln(b+a))
            = 0.5*(X+P)

            == ln(b-a)-0.5ln(b-a)-0.5ln(b+a)
            == 0.5ln(b-a)-0.5ln(b+a)
            == 0.5*(ln(b-a)-ln(b+a))

            cosh(-i*atan2(ia,b)) = b/sqrt(b^2 - a^2)
            ln(b/sqrt(b^2 - a^2)) = ln(b)-0.5*ln(b^2 - a^2)
            = ln(b) - 0.5(ln(b-a) + ln(b+a))
            = ln(b) - 0.5(X+P) = ln(cosh(z))


            ln(b+aj) = ln(b)-ln(cosh(z))+jz
            0.5((X+P)-j(X-P))
        */
        Tessarine a = t.right;
        Tessarine b = t.left;
        /*if (b.IsZero()) // Handle the special case where the left half is all zeroes (This creates problems annoyingly)
        {
            // You might have noticed that LogRecursive takes in two arguments.
            // The second argument is literally just for debugging purposes
            // So that I can track where it goes bad if it goes bad
            return LogRecursive(a, 1304) + new Vect[] { new Vect(halfpi, 1), new Vect(-halfpi, 1 + (t.dimensions >> 1)) } + shiftLogResult;
        }*/


        /*Tessarine bMa = b - a,
                  bPa = b + a;
		
		// ### This whole section is commented out because it turns out to be an attempted solution to a mathematically impossible problem.


        if ((bMa.left == bMa.right || bMa.left == -bMa.right || bPa.left == bPa.right || bPa.left == -bPa.right))
        {
			if (debug != -31415)
			{
				Vect[] vct = t.vectors;
				int indexBiggest = 0;
				double maxMag = 0;
				for (int i=0; i<vct.Length; i++)
				{
					if (Math.Abs(vct[i].value) > maxMag)
					{
						maxMag = Math.Abs(vct[i].value);
						indexBiggest = i;
					}
				}
				Vect vct1 = vct[indexBiggest];
				vct1.value = BitConverter.Int64BitsToDouble(1+BitConverter.DoubleToInt64Bits(vct1.value));
				vct[indexBiggest] = vct1;
				return LogRecursive(new Tessarine(vct), -31415);
			}
            // Basically if either b-a or b+a would be an impossible (A number that you can't divide by or take the ln of),
            // use the old-school method to get around that issue.
			// Oh also by the way, this shit doesn't work; Apparently this is impossible to fix. So there's that.
            Tessarine ath = atanh(a / b);
            return LogRecursive(b / cosh(ath)) + new Vect(1, t.dimensions >> 1) * ath + shiftLogResult;
        }*/


        Tessarine /*lna, lnb,*/ X, P;
        //lna = LogRecursive(a);
        //lnb = LogRecursive(b);
        X = LogRecursive(b - a);
        P = LogRecursive(b + a);
        // I was able to reduce every ln in the recursive definition down to just these two
        // Using the other recursive definitions, I kept getting approaching near-impossibles for high-dimensional tessarines,
        // And for doubles, those eventually round to actual impossibles
        return 0.5 * ((X + P) + new Vect(-1, t.dimensions >> 1) * (X - P)) + shiftLogResult;
        // The magic recursive formula that took me multiple days, give or take some days, to derive
    }

    // Log base 'b' of 'a', in that order
    public static Tessarine log(Tessarine b, Tessarine a)
    {
        return ln(a) / ln(b);
    }
    public static HybridTessarine log(HybridTessarine b, HybridTessarine a)
    {
        if (b.genealogy != a.genealogy) throw new InvalidOperationException("Cannot take the natural log of differing genealogies");
        return new HybridTessarine(b.genealogy, log((Tessarine)b, (Tessarine)a));
    }

    // Raise 'a' to the 'b'
    public static Tessarine pow(Tessarine a, Tessarine b)
    {
        return exp(b * ln(a));
    }
    public static HybridTessarine pow(HybridTessarine a, HybridTessarine b)
    {
        if (b.genealogy != a.genealogy) throw new InvalidOperationException("Cannot exponentiate between differing genealogies");
        return new HybridTessarine(b.genealogy, log((Tessarine)a, (Tessarine)b));
    }

    // Returns [e^t, -e^t] with some performance optimization
    public static Tessarine[] ExpNegExpPair(Tessarine t)
    {
        Tessarine a = new Tessarine(1);
        Tessarine b = new Tessarine(1);

        Vect[] v = t.vectors;
        for (int i = 0; i < v.Length; i++)
        {
            Vect[] ex = expSplit(v[i]);
            a *= ex;
            b *= new Vect[] { ex[0], -ex[1] };
        }
        return new Tessarine[] { a, b };
    }

    // Returns the square-root of 't'. Has some issues, maybe I'll fix them eventually.
    // You'll only really run into those issues if you're *trying* to
    public static Tessarine sqrt(Tessarine t)
    {
        if (t.left == t.right && t.dimensions > 2)
        {
            Tessarine sq = sqrt(0.5 * t.left);
            return sq + new Vect(1, t.dimensions >> 1) * sq;
        }
        if (t.left == -t.right && t.dimensions > 2)
        {
            Tessarine sq = sqrt(0.5 * t.left);
            return sq + new Vect(-1, t.dimensions >> 1) * sq;
        }
        return exp(0.5 * ln(t));
    }
    public static HybridTessarine sqrt(HybridTessarine t)
    {
        return new HybridTessarine(t.genealogy, sqrt((Tessarine)t));
    }

    // 2-argument arctan
    public static Tessarine atan2(Tessarine a, Tessarine b)
    {
        return -ln((b + a.imul()) * exp(-0.5 * ln(a * a + b * b))).imul();
    }
    public static HybridTessarine atan2(HybridTessarine a, HybridTessarine b)
    {
        if (b.genealogy != a.genealogy) throw new InvalidOperationException("Cannot perform atan2 on two tessarines of differing genealogies");
        return new HybridTessarine(b.genealogy, atan2((Tessarine)a, (Tessarine)b));
    }

    // Unused implementation of atan2 that was going to be used for ln,
    // But it kept giving me issues so I dropped it
    private static Tessarine[] atanh2_INTERNAL(Tessarine a, Tessarine b)
    {
        return new Tessarine[] { LogRecursive(b - a, 4), LogRecursive(b + a, 6) };
    }

    // Divide 'a' by 'b' ; a/b
    public static Tessarine div(Tessarine a, Tessarine b)
    {
        if (b.IsZeroDivisor())
        {
            return Tessarine.NaN;
        }
        return a * exp(-ln(b));
    }
    public static HybridTessarine div(HybridTessarine a, HybridTessarine b)
    {
        if (b.genealogy != a.genealogy) throw new InvalidOperationException("Cannot perform division between tessarines of differing genealogies");
        return new HybridTessarine(b.genealogy, div((Tessarine)a, (Tessarine)b));
    }

    // Honestly I don't even remember when I planned to use this,
    // Other than it was for ln. Debugging purposes, ig.
    private static Tessarine divB(Tessarine a, Tessarine b)
    {
        return a * exp(-LogRecursive(b, 5));
    }

    // Compute the cosine of Tessarine 't'
    public static Tessarine cos(Tessarine t)
    {
        Tessarine[] ex = ExpNegExpPair(t.imul());
        return 0.5 * (ex[0] + ex[1]);
    }
    public static HybridTessarine cos(HybridTessarine t)
    {
        return new HybridTessarine(t.genealogy, cos((Tessarine)t));
    }

    // Compute the sine of Tessarine 't'.
    public static Tessarine sin(Tessarine t)
    {
        Tessarine[] ex = ExpNegExpPair(t.imul());
        return -0.5 * (ex[0] - ex[1]).imul();
    }
    public static HybridTessarine sin(HybridTessarine t)
    {
        return new HybridTessarine(t.genealogy, sin((Tessarine)t));
    }

    // Compute the tangent of Tessarine 't'.
    public static Tessarine tan(Tessarine t)
    {
        Tessarine[] ex = ExpNegExpPair(t.imul());
        return -((ex[0] - ex[1]) / (ex[0] + ex[1])).imul();
    }
    public static HybridTessarine tan(HybridTessarine t)
    {
        return new HybridTessarine(t.genealogy, tan((Tessarine)t));
    }

    // Compute the cotangent of Tessarine 't'.
    public static Tessarine cot(Tessarine t)
    {
        Tessarine[] ex = ExpNegExpPair(t.imul());
        return ((ex[0] + ex[1]) / (ex[0] - ex[1])).imul();
    }
    public static HybridTessarine cot(HybridTessarine t)
    {
        return new HybridTessarine(t.genealogy, cot((Tessarine)t));
    }

    // Compute the secant of Tessarine 't'.
    public static Tessarine sec(Tessarine t)
    {
        return 1 / cos(t);
    }
    public static HybridTessarine sec(HybridTessarine t)
    {
        return new HybridTessarine(t.genealogy, sec((Tessarine)t));
    }

    // Compute the cosecant of Tessarine 't'.
    public static Tessarine csc(Tessarine t)
    {
        return 1 / sin(t);
    }
    public static HybridTessarine csc(HybridTessarine t)
    {
        return new HybridTessarine(t.genealogy, csc((Tessarine)t));
    }

    // Compute the inverse-cosine of Tessarine 't'.
    public static Tessarine acos(Tessarine t)
    {
        return halfpi + (ln(t.imul() + sqrt(1 - t * t))).imul();
    }
    public static HybridTessarine acos(HybridTessarine t)
    {
        return new HybridTessarine(t.genealogy, acos((Tessarine)t));
    }

    // Compute the inverse-sine of Tessarine 't'.
    public static Tessarine asin(Tessarine t)
    {
        return -(ln(t.imul() + sqrt(1 - t * t))).imul();
    }
    public static HybridTessarine asin(HybridTessarine t)
    {
        return new HybridTessarine(t.genealogy, asin((Tessarine)t));
    }

    // Compute the inverse-tangent of Tessarine 't'.
    public static Tessarine atan(Tessarine t)
    {
        return 0.5 * (ln((1 - t.imul())) - ln((1 + t.imul()))).imul();
    }
    public static HybridTessarine atan(HybridTessarine t)
    {
        return new HybridTessarine(t.genealogy, atan((Tessarine)t));
    }

    // Compute the inverse-cotangent of Tessarine 't'.
    public static Tessarine acot(Tessarine t)
    {
        Tessarine r1 = (1 / t).imul();
        return 0.5 * (ln((1 - r1)) - ln((1 + r1))).imul();
    }
    public static HybridTessarine acot(HybridTessarine t)
    {
        return new HybridTessarine(t.genealogy, acot((Tessarine)t));
    }

    // Compute the inverse-secant of Tessarine 't'.
    public static Tessarine asec(Tessarine t)
    {
        Tessarine r1 = 1 / t;
        return halfpi + (ln(r1.imul() + sqrt(1 - r1 * r1))).imul();
    }
    public static HybridTessarine asec(HybridTessarine t)
    {
        return new HybridTessarine(t.genealogy, asec((Tessarine)t));
    }

    // Compute the inverse-cosecant of Tessarine 't'.
    public static Tessarine acsc(Tessarine t)
    {
        Tessarine r1 = 1 / t;
        return -(ln(r1.imul() + sqrt(1 - r1 * r1))).imul();
    }
    public static HybridTessarine acsc(HybridTessarine t)
    {
        return new HybridTessarine(t.genealogy, acsc((Tessarine)t));
    }

    // Compute the hyperbolic-cosine of Tessarine 't'.
    public static Tessarine cosh(Tessarine t)
    {
        Tessarine[] ex = ExpNegExpPair(t);
        return 0.5 * (ex[0] + ex[1]);
    }
    public static HybridTessarine cosh(HybridTessarine t)
    {
        return new HybridTessarine(t.genealogy, cosh((Tessarine)t));
    }

    // Compute the hyperbolic-sine of Tessarine 't'.
    public static Tessarine sinh(Tessarine t)
    {
        Tessarine[] ex = ExpNegExpPair(t);
        return 0.5 * (ex[0] - ex[1]);
    }
    public static HybridTessarine sinh(HybridTessarine t)
    {
        return new HybridTessarine(t.genealogy, sinh((Tessarine)t));
    }

    // Compute the hyperbolic-tangent of Tessarine 't'.
    public static Tessarine tanh(Tessarine t)
    {
        Tessarine[] ex = ExpNegExpPair(t);
        return ((ex[0] - ex[1]) / (ex[0] + ex[1]));
    }
    public static HybridTessarine tanh(HybridTessarine t)
    {
        return new HybridTessarine(t.genealogy, tanh((Tessarine)t));
    }

    // Compute the hyperbolic-cotangent of Tessarine 't'.
    public static Tessarine coth(Tessarine t)
    {
        Tessarine[] ex = ExpNegExpPair(t);
        return ((ex[0] + ex[1]) / (ex[0] - ex[1]));
    }
    public static HybridTessarine coth(HybridTessarine t)
    {
        return new HybridTessarine(t.genealogy, coth((Tessarine)t));
    }

    // Compute the hyperbolic-secant of Tessarine 't'.
    public static Tessarine sech(Tessarine t)
    {
        return 1 / cosh(t);
    }
    public static HybridTessarine sech(HybridTessarine t)
    {
        return new HybridTessarine(t.genealogy, sech((Tessarine)t));
    }

    // Compute the hyperbolic-cosecant of Tessarine 't'.
    public static Tessarine csch(Tessarine t)
    {
        return 1 / sinh(t);
    }
    public static HybridTessarine csch(HybridTessarine t)
    {
        return new HybridTessarine(t.genealogy, csch((Tessarine)t));
    }

    // Compute the inverse hyperbolic-cosine of Tessarine 't'.
    public static Tessarine acosh(Tessarine t)
    {
        return ln(t + sqrt(t * t - 1));
    }
    public static HybridTessarine acosh(HybridTessarine t)
    {
        return new HybridTessarine(t.genealogy, acosh((Tessarine)t));
    }

    // Compute the inverse hyperbolic-sine of Tessarine 't'.
    public static Tessarine asinh(Tessarine t)
    {
        return ln(t + sqrt(t * t + 1));
    }
    public static HybridTessarine asinh(HybridTessarine t)
    {
        return new HybridTessarine(t.genealogy, asinh((Tessarine)t));
    }

    // Compute the inverse hyperbolic-tangent of Tessarine 't'.
    public static Tessarine atanh(Tessarine t)
    {
        return 0.5 * (ln((1 + t)) - ln((1 - t)));
    }
    public static HybridTessarine atanh(HybridTessarine t)
    {
        return new HybridTessarine(t.genealogy, atanh((Tessarine)t));
    }

    // Compute the inverse hyperbolic-cotangent of Tessarine 't'.
    public static Tessarine acoth(Tessarine t)
    {
        return atanh(1 / t);
    }
    public static HybridTessarine acoth(HybridTessarine t)
    {
        return new HybridTessarine(t.genealogy, acoth((Tessarine)t));
    }

    // Compute the inverse hyperbolic-secant of Tessarine 't'.
    public static Tessarine asech(Tessarine t)
    {
        return acosh(1 / t);
    }
    public static HybridTessarine asech(HybridTessarine t)
    {
        return new HybridTessarine(t.genealogy, asech((Tessarine)t));
    }

    // Compute the hyperbolic-cosecant of Tessarine 't'.
    public static Tessarine acsch(Tessarine t)
    {
        return asinh(1 / t);
    }
    public static HybridTessarine acsch(HybridTessarine t)
    {
        return new HybridTessarine(t.genealogy, acsch((Tessarine)t));
    }

    // Compute the conjugate of Tessarine 't'. Fuck the papers that say it doesn't exist.
    public static Tessarine conj(Tessarine t)
    {
        Tessarine nl = -ln(t);
        nl[0] = -nl[0];
        return exp(nl);
    }
    public static HybridTessarine conj(HybridTessarine t)
    {
        return new HybridTessarine(t.genealogy, conj((Tessarine)t));
    }

    private const double e_neg1 = 0.3678794411714423215955237701614608674458111310317678345078368016;
    // Uses a slightly modified version of the scipy implementation of LambertW
    public static Tessarine LambertW(Tessarine t)
    {
        if (t.IsNaN())
        {
            return t;
        }
        if (t.IsZeroDivisor() && !t.IsZero())
        {
            return Tessarine.NaN;
            // I don't currently know how to deal with this, and the scipy implementation does not account for this. I'm sure it's possible, but I don't have it.
        }

        double abs = t.Abs();

        Tessarine w;

        if (abs <= e_neg1)
        {
            if (t.IsZero())
            {
                return t;
            }
            w = t;
        }
        else
        {
            if (double.IsPositiveInfinity(t[0].value))
            {
                return t;
            }
            else if (double.IsNegativeInfinity(t[0].value))
            {
                return -t;
            }
            w = ln(t);
        }

        Tessarine ew, wew, wewz, wn = Tessarine.Zero;

        for (int i = 0; i < 100; i++)
        {
            ew = exp(w);
            wew = w * ew;
            wewz = wew - t;
            wn = w - (wewz / (wew + ew - (((w + 2) * wewz) / (w + w + 2))));

            if ((wn - w).AbsMagnitude() < 0.00000001 * (wn.AbsMagnitude()))
            {
                return wn;
            }
            else
            {
                w = wn;
            }
        }
        Console.Error.WriteLine("LambertW: Failed to converge");
        return wn;
    }
    public static HybridTessarine LambertW(HybridTessarine t)
    {
        return new HybridTessarine(t.genealogy, LambertW((Tessarine)t));
    }


    private const double fibConstL = 0.4812118250596034474977589134243684231351843343856605196610181688;
    private const double fibConstRL = -0.481211825059603447497758913424368423135184334385660519661018168;
    private const double fibConstRR = 3.1415926535897932384626433832795028841971693993751058209749445923;
    private const double fibConstB = 0.4472135954999579392818347337462552470881236719223051448541794490;

    public static Tessarine Fib(Tessarine t)
    {
        return fibConstB * (exp(t * fibConstL) - exp(t * new Tessarine(fibConstRL, fibConstRR)));
    }
    public static HybridTessarine Fib(HybridTessarine t)
    {
        return new HybridTessarine(t.genealogy, Fib((Tessarine)t));
    }
}

// The structure for an arbitrary-dimensional Tessarine
// If I didn't fuck shit up making it
public struct Tessarine
{
    // Create a Tessarine from a list of Doubles
    public Tessarine(params double[] vects)
    {

        Vect[] arr = new Vect[vects.Length.CeilDimension()];
        for (int i = 0; i < vects.Length; i++)
        {
            arr[i] = new Vect(vects[i], i);
        }
        for (int i = vects.Length; i < arr.Length; i++)
        {
            arr[i] = new Vect(0, i);
        }
        this.vects = arr;
        this.dim = (uint)arr.Length;
    }
    public Tessarine(double[] vects, bool isProperlyFormatted)
    {
        if (isProperlyFormatted)
        {
            Vect[] arrQ = new Vect[vects.Length];
            for (int i = 0; i < arrQ.Length; i++)
            {
                arrQ[i] = new Vect(vects[i], i);
            }
            this.vects = arrQ;
            this.dim = (uint)arrQ.Length;
            return;
        }

        Vect[] arr = new Vect[vects.Length.CeilDimension()];
        for (int i = 0; i < vects.Length; i++)
        {
            arr[i] = new Vect(vects[i], i);
        }
        for (int i = vects.Length; i < arr.Length; i++)
        {
            arr[i] = new Vect(0, i);
        }
        this.vects = arr;
        this.dim = (uint)arr.Length;
    }
    // Create a Tessarine from a list of Vectors
    public Tessarine(params Vect[] vects)
    {
        this.vects = null;
        this.dim = 2;

        Vect[] arr = CreateVectArray(vects.GetHighestDimension());
        for (int i = 0; i < vects.Length; i++)
        {
            Vect v = vects[i];
            arr[v.index].value += v.value;
        }
        this.vects = arr;
        this.dim = (uint)arr.Length;
    }
    public Tessarine(Vect[] vects, bool isProperlyFormatted)
    {
        if (isProperlyFormatted)
        {
            this.vects = vects;
            this.dim = (uint)vects.Length;
            return;
        }

        this.vects = null;
        this.dim = 2;

        Vect[] arr = CreateVectArray(vects.GetHighestDimension());
        for (int i = 0; i < vects.Length; i++)
        {
            Vect v = vects[i];
            arr[v.index].value += v.value;
        }
        this.vectors = arr;
    }

    // Convert this tessarine to a String (Plain Text)
    // I have yet to implement the inverse, but I think it's better not to have it
    public override string ToString()
    {
        if (this.IsNaN())
        {
            return "NaN";
        }
        string res = this[0].ToString();
        for (int i = 1; i < this.dimensions; i++)
        {
            res += this[i].ToString(PrintFlags.ForceSign | PrintFlags.SpaceSign);
        }
        return res;
    }

    // Returns true if any value inside this Tessarine is NaN
    public bool IsNaN()
    {
        for (int i = 0; i < this.vects.Length; i++)
        {
            if (double.IsNaN(this.vects[i].value))
            {
                return true;
            }
        }
        return false;
    }

    // Get the sum of absolute-values of each vector.
    // Useful for finding the difference between two Tessarines,
    // Particularly for finding the computational error in certain operations
    public double AbsMagnitude()
    {
        double d = 0;
        for (int i = 0; i < this.vects.Length; i++)
        {
            d += Math.Abs(this.vects[i].value);
        }
        return d;
    }

    // Recursively compute absolute-value. This is necessary because without this, we run into issues with numbers being too big for large tessarines and we get infinities easily.
    private static double RecursiveAbs(double[] v, uint rot, uint len)
    {
        if (len != 1)
        {
            return Math.Sqrt(RecursiveAbs(v, rot, len >> 1) * RecursiveAbs(v, rot + len, len >> 1));
        }

        // Deepest level of the recursion tree, endpoint
        double left = 0;
        double right = 0;
        for (uint i = 0; i < v.Length; i += 2)
        {
            if (((i & rot).GetSetBitCount() & 1) == 0)
            {
                left += v[i];
                right += v[i + 1];
            }
            else
            {
                left -= v[i];
                right -= v[i + 1];
            }
        }
        left *= left;
        right *= right;
        return Math.Sqrt(left + right);
    }

    public static Tessarine Random(uint dimensions)
    {
        return Random(dimensions, 0, 1);
    }
    public static Tessarine Random(int dimensions)
    {
        return Random((uint)dimensions, 0, 1);
    }
    public static Tessarine Random(int dimensions, double min, double max)
    {
        return Random((uint)dimensions, min, max);
    }
    public static Tessarine Random(uint dimensions, double min, double max)
    {
        Vect[] v = new Vect[dimensions];
        double range = max - min;
        for (int i = 0; i < dimensions; i++)
        {
            v[i] = new Vect((Global.rand.NextDouble() * range) + min, i);
        }
        return new Tessarine(v);
    }

    // Creates a zeroed Vect array of the given length (Rounded *up* to the closest suitable dimension)
    public static Vect[] CreateVectArray(int length)
    {
        return CreateVectArray((uint)length);
    }
    // Same as above but better
    public static Vect[] CreateVectArray(uint length)
    {
        Vect[] arr = new Vect[length.CeilDimension()];
        for (int i = 0; i < arr.Length; i++)
        {
            arr[i] = new Vect(0, i);
        }
        return arr;
    }
    // Takes in a Vect array, fixes the ordering (And adds same-index Vectors together),
    // And returns a new Vect[] array with the smallest suitable dimensions
    public static Vect[] CreateVectArray(Vect[] arr)
    {
        Vect[] res = CreateVectArray(arr.GetHighestDimension());
        for (int i = 0; i < arr.Length; i++)
        {
            res[arr[i].index].value += arr[i].value;
        }
        return res;
    }
    // Resizes this Tessarine to the dimensions given (Ceiled/Rounded up to closest suitable dimension count)
    public Tessarine ToDimension(uint dimension)
    {
        dimension = dimension.CeilDimension();
        Vect[] newVects = CreateVectArray(dimension);
        uint min = dimension;
        if (this.dimensions < min) min = this.dimensions;

        for (int i = 0; i < min; i++)
        {
            newVects[i] = this[i];
        }
        this.vectors = newVects;
        return this;
    }

    // Returns true if this Tessarine is a zero-divisor.
    /*
		The criteria for a Tessarine q=a+bj being a zero-divisor is as follows:
		
		1: q is zero
		2: q is not reducable to a 2-dimensional Tessarine (Complex number) >>> If true, CONTINUE
		3: a+b is a zero-divisor (Recursive down to 2-dimensional)
		4: a-b is a zero-divisor (Recursive down to 2-dimensional)
	*/
    public bool IsZeroDivisor()
    {
        if (this.IsZero())
        {
            return true;
        }
        if (this.dimensions == 2)
        {
            return false;
        }
        Tessarine l = this.left;
        Tessarine r = this.right;
        if ((l - r).IsZeroDivisor())
        {
            return true;
        }
        if ((l + r).IsZeroDivisor())
        {
            return true;
        }
        return false;
    }

    // Reduces this Tessarine to the smallest-dimensional Tessarine that represents the exact same value.
    public Tessarine Reduce()
    {
        if (this.dimensions == 2) return this.Copy();
        Vect[] v = this.vectors;
        uint i = this.dimensions >> 1;
        for (; i < v.Length; i++)
        {
            if (v[i].value != 0)
            {
                return this.Copy();
            }
        }
        return this.left.Reduce();
    }

    // Creates a copy of this Tessarine
    public Tessarine Copy()
    {
        return new Tessarine(Tessarine.CreateVectArray(this.vectors));
    }

    // Fixes the ordering of the terms in this Tessarine if they are fucked up for some reason (You have to save the result to a variable)
    public Tessarine FixOrdering()
    {
        Vect[] arr = new Vect[this.vectors.GetHighestDimension()];
        for (int i = 0; i < arr.Length; i++)
        {
            arr[i] = new Vect(0, i);
        }
        for (int i = 0; i < this.vects.Length; i++)
        {
            Vect v = this.vects[i];
            Vect v1 = arr[v.index];
            v1.value += v.value;
            arr[v.index] = v1;
        }
        this.vectors = arr;
        return this;
    }

    // Multiply this tessarine by 'i' (First imaginary unit)
    public Tessarine imul()
    {
        Vect[] res = new Vect[this.dimensions];
        for (int i = 0; i < this.dimensions; i++)
        {
            Vect r1 = this[i].imul();
            res[r1.index] = r1;
        }
        return new Tessarine(res);
    }

    // Returns 'true' if every term in this Tessarine is zero; Otherwise, 'false'
    public bool IsZero()
    {
        for (int i = 0; i < this.vects.Length; i++)
        {
            if (this.vects[i].value != 0)
            {
                return false;
            }
        }
        return true;
    }

    // Multiply two Tessarines
    public static Tessarine operator *(Tessarine left, Tessarine right)
    {
        Vect[] mulBV = right.vectors;
        Vect[] addVals = new Vect[mulBV.Length];
        Tessarine t = new Tessarine(0, 0);
        for (int i = 0; i < left.dimensions; i++)
        {
            if (left[i].value == 0) continue;
            for (int j = 0; j < mulBV.Length; j++)
            {
                if (mulBV[j].value == 0) continue;
                addVals[j] = left[i] * mulBV[j];
            }
            t += addVals;
        }
        return t;
    }

    // Multiply a Tessarine by a Vector
    public static Tessarine operator *(Tessarine left, Vect right)
    {
        if (right.value == 0) return Tessarine.Zero;
        /*if (left.dimensions <= right.index)
        {
            left = left.ToDimension(Math.Max(left.dimensions, (right.index + 1).CeilDimension()));
        }*/
        Vect[] vectArr = left.vectors;
        Vect[] resArr = CreateVectArray((uint)(vectArr.Length - 1) | ((right.index | 1) - 1));
        for (int i = 0; i < vectArr.Length; i++)
        {
            Vect v = vectArr[i] * right;
            resArr[v.index] = v;
        }
        return new Tessarine(resArr, true);
    }

    // Multiply a Vector by a Tessarine
    public static Tessarine operator *(Vect left, Tessarine right)
    {
        return right * left;
    }

    // Multiply a Tessarine by an array of Vectors
    public static Tessarine operator *(Tessarine left, Vect[] right)
    {
        Tessarine res = new Tessarine(0);
        for (int i = 0; i < right.Length; i++)
        {
            if (right[i].value == 0) continue;
            res += left * right[i];
        }
        return res;
    }

    // Multiply a Vector array by a Tessarine
    public static Tessarine operator *(Vect[] left, Tessarine right)
    {
        return right * left;
    }

    // Add a Tessarine and a double
    public static Tessarine operator +(Tessarine left, double right)
    {
        Tessarine t = left.Copy();
        Vect v = t[0];
        v.value += right;
        t[0] = v;
        return t;
    }

    // Same as above
    public static Tessarine operator +(double left, Tessarine right)
    {
        return right + left;
    }

    // Subtract a double from a Tessarine
    public static Tessarine operator -(Tessarine left, double right)
    {
        Tessarine t = left.Copy();
        Vect v = t[0];
        v.value -= right;
        t[0] = v;
        return t;
    }

    // Subtract a Tessarine from a double
    public static Tessarine operator -(double left, Tessarine right)
    {
        Tessarine t = -right.Copy();
        Vect v = t[0];
        v.value = left + v.value;
        t[0] = v;
        return t;
    }

    // Multiply a Tessarine by a double
    public static Tessarine operator *(Tessarine left, double right)
    {
        Vect[] v = new Vect[left.dimensions];
        Vect[] vL = left.vectors;
        for (int i = 0; i < v.Length; i++)
        {
            v[i] = right * vL[i];
        }
        return new Tessarine(v, true);
    }

    // Same as above
    public static Tessarine operator *(double left, Tessarine right)
    {
        return right * left;
    }

    // Divide a Tessarine by a double;
    public static Tessarine operator /(Tessarine left, double right)
    {
        Vect[] v = left.vectors.Copy();
        right = 1 / right;
        for (int i = 0; i < v.Length; i++)
        {
            v[i] *= right;
        }
        return new Tessarine(v);
    }

    // Divide a double by a Tessarine
    public static Tessarine operator /(double left, Tessarine right)
    {
        if (left == 1)
        {
            return TessMath.exp(-TessMath.ln(right));
        }
        else if (left == -1)
        {
            return -TessMath.exp(-TessMath.ln(right));
        }
        return left * TessMath.exp(-TessMath.ln(right));
    }

    // Divide a Vector by a Tessarine
    public static Tessarine operator /(Vect left, Tessarine right)
    {
        return left * TessMath.exp(-TessMath.ln(right));
    }

    // Divide a Tessarine by a Vector
    public static Tessarine operator /(Tessarine left, Vect right)
    {
        Vect[] v = left.vectors;
        for (int i = 0; i < v.Length; i++)
        {
            v[i] = v[i] / right;
        }
        return new Tessarine(v);
    }

    // Divide a Tessarine by a Tessarine
    public static Tessarine operator /(Tessarine left, Tessarine right)
    {
        return TessMath.div(left, right);
    }

    // DONE # Todo: Add ln operation
    // DONE # Todo: Add conj operation
    // DONE # Todo: Add division operation
    // DONE # Todo: Add sqrt operation
    // Update: Not feasible --- Todo: Add Decimal version of Tessarines for higher precision
    // Todo: Add LambertW
    // Todo: Make Random Tessarine easier to make


    // Add a Tessarine and a Vector array
    public static Tessarine operator +(Tessarine left, Vect[] right)
    {
        //Tessarine t = left.Copy();
        //t = t.ToDimension(Math.Max(left.dimensions, right.GetHighestDimension()));
        Vect[] tVects = left.vectors;
        Vect[] res = new Vect[((left.dimensions - 1) | (right.GetHighestDimension() - 1)).CeilDimension()];
        for (int i = 0; i < tVects.Length; i++)
        {
            res[i] = tVects[i];
        }
        for (int i = tVects.Length; i < res.Length; i++)
        {
            res[i] = new Vect(0, i);
        }
        for (int i = 0; i < right.Length; i++)
        {
            Vect r = right[i];
            res[r.index].value += r.value;
        }
        return new Tessarine(res, true);
        /*for (int i = 0; i < right.Length; i++)
        {
            Vect r = right[i];
            Vect v = t[r.index];
            v.value += r.value;
            t[r.index] = v;
        }
        return t;*/
    }

    // Same as above
    public static Tessarine operator +(Vect[] left, Tessarine right)
    {
        return right + left;
    }

    // Add a Tessarine to another Tessarine
    public static Tessarine operator +(Tessarine left, Tessarine right)
    {
        return left + right.vectors;
    }

    // Just copies the Tessarine value
    public static Tessarine operator +(Tessarine val)
    {
        return val.Copy();
    }

    // Subtract a Vector array from a Tessarine
    public static Tessarine operator -(Tessarine left, Vect[] right)
    {
        Tessarine t = left.Copy();
        t = t.ToDimension(Math.Max(left.dimensions, right.GetHighestDimension()));

        for (int i = 0; i < right.Length; i++)
        {
            Vect r = right[i];
            Vect v = t[r.index];
            v.value -= r.value;
            t[r.index] = v;
        }
        return t;
    }

    // Subtract a Tessarine from a Vector array
    public static Tessarine operator -(Vect[] left, Tessarine right)
    {
        Tessarine t = right.Copy();
        t = t.ToDimension(Math.Max(right.dimensions, left.GetHighestDimension()));

        for (int i = 0; i < left.Length; i++)
        {
            Vect r = left[i];
            Vect v = t[r.index];
            v.value = r.value - v.value;
            t[r.index] = v;
        }
        return t;
    }

    // Subtract a Tessarine from a Tessarine
    public static Tessarine operator -(Tessarine left, Tessarine right)
    {
        return left - right.vectors;
    }

    // Negate this Tessarine
    public static Tessarine operator -(Tessarine val)
    {
        Vect[] vects = val.vectors.Copy();
        for (int i = 0; i < vects.Length; i++)
        {
            vects[i] = -vects[i];
        }
        return new Tessarine(vects);
    }

    // Exponentiate a Tessarine by another Tessarine
    public static Tessarine operator ^(Tessarine left, Tessarine right)
    {
        return TessMath.pow(left, right);
    }

    // Take the sqrt of this Tessarine
    public Tessarine sqrt()
    {
        return TessMath.sqrt(this);
    }

    // Square this Tessarine
    public Tessarine sq()
    {
        return this * this;
    }

    // Get the reciprocal of this Tessarine
    public Tessarine reciprocal()
    {
        return 1 / this;
    }

    // Return the exp (e^x) of this Tessarine
    public Tessarine exp()
    {
        return TessMath.exp(this);
    }

    // Compute the natural log of this Tessarine
    public Tessarine ln()
    {
        return TessMath.ln(this);
    }

    // Returns the absolute-value of this Tessarine
    public double Abs()
    {
        Vect[] v = this.FixOrdering().vectors;
        uint dimensions = this.dimensions;

        double[] vs = new double[v.Length];
        for (int i = 0; i < vs.Length; i++)
        {
            vs[i] = v[i].value;
        }

        return RecursiveAbs(vs, 0, dimensions >> 1);
    }

    public Tessarine conj()
    {
        return TessMath.conj(this);
    }

    // Check for value-equality between both Tessarines
    public static bool operator ==(Tessarine left, Tessarine right)
    {
        left = left.Reduce();
        right = right.Reduce();
        if (left.dimensions != right.dimensions)
        {
            return false;
        }
        for (int i = 0; i < left.vects.Length; i++)
        {
            if (left[i] != right[i])
            {
                return false;
            }
        }
        return true;
    }

    // Check for value-inequality between both Tessarines
    public static bool operator !=(Tessarine left, Tessarine right)
    {
        return !(left == right);
    }

    // Get or set the vector at the given index
    // This is where you can fuck up the ordering on the 'set' part
    public Vect this[int index]
    {
        get
        {
            if (index >= this.vectors.Length)
            {
                if (this.IsNaN())
                {
                    return new Vect(double.NaN, index);
                }
                return new Vect(0, index);
            }
            return this.vectors[index];
        }
        set
        {
            if (index >= this.vectors.Length)
            {
                Vect[] arr = CreateVectArray(index);
                for (int i = 0; i < this.vectors.Length; i++)
                {
                    Vect v = arr[this.vectors[i].index];
                    v.value += this.vectors[i].value;
                    arr[this.vectors[i].index] = v;
                }
                this.vectors = arr;
            }
            this.vectors[index] = value;
        }
    }

    // Same as above but no negative indexes
    public Vect this[uint index]
    {
        get
        {
            if (index >= this.vectors.Length)
            {
                if (this.IsNaN())
                {
                    return new Vect(double.NaN, index);
                }
                return new Vect(0, index);
            }
            return this.vectors[(int)index];
        }
        set
        {
            if (index >= this.vectors.Length)
            {
                Vect[] arr = CreateVectArray(index);
                for (int i = 0; i < this.vectors.Length; i++)
                {
                    Vect v = arr[this.vectors[i].index];
                    v.value += this.vectors[i].value;
                    arr[this.vectors[i].index] = v;
                }
                this.vectors = arr;
            }
            this.vectors[(int)index] = value;
        }
    }

    // Implicitly convert a Vector to a Tessarine
    public static implicit operator Tessarine(Vect v)
    {
        return new Tessarine(v);
    }

    // Implicitly convert a Double to a Tessarine
    public static implicit operator Tessarine(double d)
    {
        return new Tessarine(d);
    }

    // Implicitly convert a Vector array to a Tessarine
    public static implicit operator Tessarine(Vect[] v)
    {
        return new Tessarine(v);
    }

    // Implicitly convert a Double array to a Tessarine
    public static implicit operator Tessarine(double[] d)
    {
        return new Tessarine(d);
    }

    // Is this Tessarine nonzero?
    public static bool operator true(Tessarine t)
    {
        return !t.IsZero();
    }

    // Is this Tessarine zero?
    public static bool operator false(Tessarine t)
    {
        return t.IsZero();
    }

    // Implicitly convert a Tessarine to a Vector array
    public static implicit operator Vect[](Tessarine t)
    {
        return t.vectors.Copy();
    }

    // Implicitly convert a Tessarine to a Double array
    public static implicit operator double[](Tessarine t)
    {
        double[] d = new double[t.vects.Length];
        for (int i = 0; i < d.Length; i++)
        {
            d[i] = t.vects[i].value;
        }
        return d;
    }

    // Get the left half of this Tessarine
    public Tessarine left
    {
        get
        {
            Vect[] v = new Vect[this.dimensions >> 1];
            for (int i = 0; i < v.Length; i++)
            {
                v[i] = this.vectors[i];
            }
            return new Tessarine(v);
        }
    }

    // Get the right half of this Tessarine
    public Tessarine right
    {
        get
        {
            uint p = this.dimensions >> 1;
            Vect[] v = new Vect[p];
            for (int i = 0; i < v.Length; i++)
            {
                Vect v1 = this.vectors[v.Length + i];
                v1.index ^= p;
                v[i] = v1;
            }
            return new Tessarine(v);
        }
    }

    // Get or set the vectors in this Tessarine
    public Vect[] vectors
    {
        get
        {
            return this.vects;
        }
        set
        {
            this.vects = CreateVectArray(value);
            this.dimensions = (uint)this.vects.Length;
        }
    }

    // Get the dimensions of this Tessarine
    public uint dimensions
    {
        get
        {
            return dim;
        }
        private set
        {
            dim = value;
        }
    }
    private Vect[] vects;
    private uint dim;

    // Get the default Zero tessarine [0,0]
    public static Tessarine Zero
    {
        get
        {
            return new Tessarine(0, 0);
        }
    }

    // Returns a NaN valued Tessarine
    public static Tessarine NaN
    {
        get
        {
            return new Tessarine(double.NaN, double.NaN);
        }
    }
}

// Class containing extension methods to make life easier
public static class Global
{
    public static Random rand = new Random();
    // Returns the number of set bits in the binary representation of this unsigned integer
    public static int GetSetBitCount(this uint val)
    {
        uint res = 0;
        while (val != 0)
        {
            res += val & 1;
            val >>= 1;
        }
        return (int)res;
    }
    // Returns the number of set bits in the binary representation of this integer
    public static int GetSetBitCount(this int val)
    {
        return GetSetBitCount((uint)val);
    }

    // Gets the index of the highest set bit
    public static int HighestSetBitIndex(this int val)
    {
        int res = 0;
        uint val0 = ((uint)val) >> 1;
        while (val0 > 0)
        {
            res++;
            val >>= 1;
        }
        return res;
    }
    public static int HighestSetBitIndex(this uint val)
    {
        int res = 0;
        uint val0 = (val) >> 1;
        while (val0 > 0)
        {
            res++;
            val >>= 1;
        }
        return res;
    }
    public static int LowestSetBitIndex(this int val)
    {
        int res = 0;
        uint val0 = (uint)val;
        while (val != 0)
        {
            if ((val0 & 1) == 1)
            {
                return res;
            }
            res++;
            val0 >>= 1;
        }
        return 32;
    }
    public static int LowestSetBitIndex(this uint val)
    {
        int res = 0;
        uint val0 = val;
        while (val != 0)
        {
            if ((val0 & 1) == 1)
            {
                return res;
            }
            res++;
            val0 >>= 1;
        }
        return 32;
    }

    public static int SwapBits(this int val, int pos1, int pos2)
    {
        uint uval = (uint)val;
        if (pos1 < pos2)
        {
            pos1 ^= pos2;
            pos2 ^= pos1;
            pos1 ^= pos2;
        }
        uint mask1 = 1u << pos1;
        uint mask2 = 1u << pos2;
        uint val1 = uval & mask1;
        uint val2 = uval & mask2;

        uint mask = 4294967295 - (mask1 | mask2);

        int diff = pos1 - pos2;

        uval = (uval & mask) | (val1 >> diff) | (val2 << diff);
        return (int)uval;
    }
    public static uint SwapBits(this uint uval, int pos1, int pos2)
    {
        if (pos1 < pos2)
        {
            pos1 ^= pos2;
            pos2 ^= pos1;
            pos1 ^= pos2;
        }
        uint mask1 = 1u << pos1;
        uint mask2 = 1u << pos2;
        uint val1 = uval & mask1;
        uint val2 = uval & mask2;

        uint mask = 4294967295 - (mask1 | mask2);

        int diff = pos1 - pos2;

        uval = (uval & mask) | (val1 >> diff) | (val2 << diff);
        return uval;
    }

    // Rounds this int to the closest 2^n equal to or higher than the current value
    public static int CeilDimension(this int val)
    {
        int res = 1;
        uint val0 = ((uint)val - 1) >> 1;
        while (val0 > 0)
        {
            res++;
            val0 >>= 1;
        }
        return 1 << (res);
    }
    // Same as above
    public static uint CeilDimension(this uint val)
    {
        uint res = 1;
        uint val0 = (val - 1) >> 1;
        while (val0 > 0)
        {
            res++;
            val0 >>= 1;
        }
        return (uint)1 << (int)(res);
    }

    public static bool IsAllSameGenealogy(this HybridVect[] val)
    {
        for (int i = 0; i < val.Length - 1; i++)
        {
            if (val[i].genealogy != val[i + 1].genealogy) return false;
        }
        return true;
    }

    // Gets the highest dimension of the Vectors in this array (Rounded up to 2^n)
    public static uint GetHighestDimension(this Vect[] val)
    {
        uint max = 0;
        for (int i = 0; i < val.Length; i++)
        {
            uint c = val[i].index;
            if (c > max)
            {
                max = c;
            }
        }
        return (max + 1).CeilDimension();
    }
    public static uint GetHighestDimension(this HybridVect[] val)
    {
        uint max = 0;
        for (int i = 0; i < val.Length; i++)
        {
            uint c = val[i].index;
            if (c > max)
            {
                max = c;
            }
        }
        return (max + 1).CeilDimension();
    }

    // Creates a copy of this Vector array
    public static Vect[] Copy(this Vect[] arr)
    {
        return Tessarine.CreateVectArray(arr);
    }
    public static HybridVect[] Copy(this HybridVect[] arr)
    {
        return HybridTessarine.CreateVectArray(arr);
    }

    // Returns a string of this Double value rounded to 11 significant figures
    public static string RoundString(this double d)
    {
        return d.RoundString(11);
    }

    // Returns a string of this Double value rounded to a number of significant figures of your choosing
    public static string RoundString(this double d, int digits)
    {
        string dt = d.ToString();
        string e = "";
        string sig = "";

        int indE = dt.IndexOf("E");
        if (indE != -1)
        {
            e = dt.Substring(indE);
            sig = dt.Substring(0, indE);
        }
        else
        {
            sig = dt;
        }
        if (sig.Length < digits + 1) return dt;
        int ind = sig.IndexOf(".");
        if (ind == -1) return dt;
        if (ind > digits) return ((long)(Math.Round(d))).ToString();
        sig = sig.Replace(".", "");
        string mSign = "";
        int countAppendZeroes = 0;
        if (sig.StartsWith("-"))
        {
            mSign = "-";
            sig = sig.Substring(1);
        }
        int startLen = sig.Length;
        sig = sig.TrimStart('0');
        countAppendZeroes = startLen - sig.Length;
        if (sig.Length < digits) return dt;
        sig = sig.Insert(digits, ".");
        sig = Math.Round(double.Parse(sig)).ToString();
        if (sig.Length > digits)
        {
            ind++;
            countAppendZeroes--;
        }
        if (countAppendZeroes > 0)
        {
            sig = new String('0', countAppendZeroes) + sig;
        }
        sig = mSign + sig;
        if (sig.Length > ind)
        {
            sig = sig.Replace(".", "").Insert(ind, ".").TrimEnd('0');
        }
        else
        {
            sig = sig.Replace(".", "");
        }
        if (sig.EndsWith("."))
        {
            sig = sig.TrimEnd('.');
        }
        return sig + e;
    }
}

// Flags for how to print Vectors to strings. Some are unused.
// I'll add this to Tessarines at some point.
[Flags]
public enum PrintFlags
{
    None = 0,
    HideZeroes = 1,
    Reorder = 2,
    ForceSign = 4,
    MatrixForm = 8,
    SpaceSign = 16
}

// A single Vector, with a magnitude and a Vector index
public struct Vect
{
    public double value; // The magnitude of this Vector (Can be negative)
    public uint index; // The index of this vector

    public static readonly Vect i = new Vect(1, 1); // The first imaginary vector
    public static int DefaultSigfigs = 9;
    // The number of significant figures to round to when printing this Vector to a String. This can be modified at runtime with code.

    // The subscript values that are used to indicate the Vector index.
    // This particular string wasn't utilized in the code, but it could've been.
    private const string SubscriptValues = "₀₁₂₃₄₅₆₇₈₉";

    // Convert this Vector to a String.
    public override string ToString()
    {
        string s = this.value.RoundString(DefaultSigfigs);
        if (s.Contains("E"))
        {
            s = "(" + s + ")";
        }
        return s + "e" + ToSubscript(this.index);
    }
    // Convert this Vector to a String, with the provided rules.
    // PrintFlags.MatrixForm and PrintFlags.Reorder are currently not supported,
    // But PrintFlags.ForceSign, PrintFlags.HideZeroes, and PrintFlags.SpaceSign are
    public string ToString(PrintFlags flags)
    {
        string res = "";
        if ((flags & PrintFlags.MatrixForm) != 0)
        {
            throw new NotImplementedException();
        }
        if ((flags & PrintFlags.Reorder) != 0)
        {
            throw new NotImplementedException();
        }
        if ((flags & PrintFlags.HideZeroes) != 0 && this.value == 0)
        {
            return "";
        }
        bool isNegative = this.value < 0;
        string s = (isNegative ? -this.value : this.value).RoundString(DefaultSigfigs);
        if (s.Contains("E"))
        {
            s = "(" + s + ")";
        }
        if ((flags & PrintFlags.ForceSign) != 0 && (this.value >= 0 || s.Contains("E")))
        {
            res += "+";
        }
        else if (this.value < 0)
        {
            res += "-";
        }

        if ((flags & PrintFlags.SpaceSign) != 0)
        {
            res = " " + res + " ";
        }
        res += s + "e" + ToSubscript(this.index);
        return res;
    }

    // Convert a uint to a subscript string.
    private static string ToSubscript(uint num)
    {
        if (num == 0)
        {
            return "₀";
        }
        string res = "";
        while (num > 0)
        {
            char val = (char)(num % 10);
            num /= 10;
            res = ((char)(val + '₀')).ToString() + res;
        }
        return res;
    }

    // Create a Vector from a magnitude and an index
    public Vect(double val, uint ind)
    {
        this.value = val;
        this.index = ind;
    }

    // Create a Vector from a magnitude and an index
    // Be advised: -1 in here is index 4294967295,
    // Which will almost certainly crash when you try to make a Tessarine out of that
    public Vect(double val, int ind)
    {
        this.value = val;
        this.index = (uint)ind;
    }

    // Multiply this Vector by 'i', the first imaginary unit
    public Vect imul()
    {
        if ((this.index & 1) == 0)
        {
            return new Vect(-this.value, this.index ^ 1);
        }
        return new Vect(this.value, this.index ^ 1);
    }

    // Does nothing but prevent an error
    public static Vect operator +(Vect v)
    {
        return v;
    }

    // Negate this Vector
    public static Vect operator -(Vect v)
    {
        v.value = -v.value;
        return v;
    }

    // Multiply two Vectors
    public static Vect operator *(Vect left, Vect right)
    {
        double co = left.value * right.value;
        uint ind = left.index ^ right.index;
        if (((left.index & right.index) & 1) == 0)
        {
            return new Vect(co, ind);
        }
        return new Vect(-co, ind);
    }

    // Multiply a Vector and a double
    public static Vect operator *(double left, Vect right)
    {
        right.value *= left;
        return right;
    }
    // Same as above
    public static Vect operator *(Vect left, double right)
    {
        left.value *= right;
        return left;
    }

    // Divide a Vector by a double
    public static Vect operator /(Vect left, double right)
    {
        left.value /= right;
        return left;
    }

    // Divide a double by a Vector
    public static Vect operator /(double left, Vect right)
    {
        double co = left / right.value;
        if ((right.index & 1) == 0)
        {
            return new Vect(co, right.index);
        }
        return new Vect(-co, right.index);
    }

    // Divide a Vector by another Vector
    public static Vect operator /(Vect left, Vect right)
    {
        double resD = left.value / right.value;
        if ((right.index & 1) == 1 && (left.index & 1) == 0)
        {
            return new Vect(-resD, left.index ^ right.index);
        }
        return new Vect(resD, left.index ^ right.index);
    }

    // Add two Vectors together
    public static Tessarine operator +(Vect left, Vect right)
    {
        return new Tessarine(left, right);
    }

    // Subtract two Vectors
    public static Tessarine operator -(Vect left, Vect right)
    {
        return new Tessarine(left, -right);
    }

    // Implicitly convert a Vector to a double (Loses index information)
    public static implicit operator double(Vect vec)
    {
        return vec.value;
    }

    // Implicitly convert a Vector to an integer (Returns index, loses magnitude information)
    public static explicit operator uint(Vect vec)
    {
        return vec.index;
    }

    // Is this vector nonzero
    public static bool operator true(Vect vec)
    {
        return vec.value != 0;
    }

    // Is this vector zero
    public static bool operator false(Vect vec)
    {
        return vec.value == 0;
    }

    // Compares two Vectors for equality
    public static bool operator ==(Vect left, Vect right)
    {
        return (left.index == right.index) && (left.value == right.value);
    }

    // Compares two Vectors for inequality
    public static bool operator !=(Vect left, Vect right)
    {
        return !(left == right);
    }
}

public static class Genealogy
{
    public static uint Generate(params int[] genealogy_flags)
    {
        uint res = 0;
        for (int i = 0; i < genealogy_flags.Length; i++)
        {
            if (genealogy_flags[i] != -1 && genealogy_flags[i] != +1) throw new ArgumentException("The provided genealogies for Genealogy.Generate() must be either -1 for rotators, or +1 for reflectors.");
            res |= (uint)(genealogy_flags[i] == -1 ? 1 : 0) << i;
        }
        return res;
    }
}

public struct HybridVect
{
    public uint genealogy;
    public uint index;
    public double value;

    public HybridVect(double value, uint index, uint genealogy)
    {
        this.value = value;
        this.index = index;
        this.genealogy = genealogy;
    }
    public HybridVect(double value, int index, int genealogy)
    {
        this.value = value;
        this.index = (uint)index;
        this.genealogy = (uint)genealogy;
    }
    public HybridVect(double value, int index, uint genealogy)
    {
        this.value = value;
        this.index = (uint)index;
        this.genealogy = genealogy;
    }
    public HybridVect(double value, uint index, int genealogy)
    {
        this.value = value;
        this.index = index;
        this.genealogy = (uint)genealogy;
    }
    public HybridVect(Vect v, uint genealogy)
    {
        if (genealogy == 0)
        {
            if ((v.index & 1) == 1)
            {
                throw new InvalidOperationException("Unable to convert a rotator to a pure reflector genealogy");
            }
            else
            {
                this.value = v.value;
                this.index = v.index >> 1;
                this.genealogy = 0;
                return;
            }
        }

        this.value = v.value;
        this.genealogy = genealogy;
        int lowestSet = genealogy.LowestSetBitIndex();
        this.index = v.index.SwapBits(lowestSet, 0);
        genealogy = genealogy.SwapBits(lowestSet, 0);

        uint bitand = (this.index & genealogy) & 4294967294;
        uint setcnt = (uint)bitand.GetSetBitCount();
        this.index ^= (setcnt & 1) << lowestSet;
        if ((((setcnt & 3) - (v.index & 1) + 1) & 2) != 0)
        {
            this.value = -this.value;
        }
    }

    public static implicit operator Vect(HybridVect val)
    {
        if (val.genealogy == 0)
        {
            val.index <<= 1;
            val.genealogy = (val.genealogy << 1) | 1;
        }

        int lowestSet = val.genealogy.LowestSetBitIndex();
        val.index = val.index.SwapBits(lowestSet, 0);
        val.genealogy = val.genealogy.SwapBits(lowestSet, 0);

        uint bitand = (val.index & val.genealogy) & 4294967294;
        uint setcnt = (uint)bitand.GetSetBitCount();
        Vect res = new Vect(val.value, val.index ^ (setcnt & 1));
        if ((((setcnt & 3) + (val.index & 1)) & 2) != 0)
        {
            res.value = -res.value;
        }
        return res;
    }
    public static implicit operator double(HybridVect val)
    {
        return val.value;
    }
    public static explicit operator uint(HybridVect val)
    {
        return val.index;
    }
    public static HybridVect operator +(HybridVect val)
    {
        return val;
    }
    public static HybridVect operator -(HybridVect val)
    {
        return new HybridVect(-val.value, val.index, val.genealogy);
    }

    public static HybridVect operator *(HybridVect left, HybridVect right)
    {
        if (left.genealogy != right.genealogy) throw new NotImplementedException("Multiplication between different genealogies is not yet supported");

        double co = left.value * right.value;
        uint ind = left.index ^ right.index;
        uint genealogy = left.genealogy;
        if ((((left.index & right.index & genealogy).GetSetBitCount()) & 1) == 0)
        {
            return new HybridVect(co, ind, genealogy);
        }
        return new HybridVect(-co, ind, genealogy);
    }
    public static HybridVect operator *(HybridVect left, double right)
    {
        left.value *= right;
        return left;
    }
    public static HybridVect operator *(double left, HybridVect right)
    {
        right.value *= left;
        return right;
    }

    public static HybridVect operator /(HybridVect left, HybridVect right)
    {
        if (left.genealogy != right.genealogy) throw new NotImplementedException("Division between different genealogies is not yet supported");
        double resD = left.value / right.value;
        uint genealogy = left.genealogy;

        if (((right.index & genealogy).GetSetBitCount() & 1) == 1 && ((left.index & right.index & genealogy).GetSetBitCount() & 1) == 0)
        {
            return new HybridVect(-resD, left.index ^ right.index, genealogy);
        }
        return new HybridVect(resD, left.index ^ right.index, genealogy);
    }
    public static HybridVect operator /(HybridVect left, double right)
    {
        left.value /= right;
        return left;
    }
    public static HybridVect operator /(double left, HybridVect right)
    {
        if (((right.index & right.genealogy).GetSetBitCount() & 1) == 1)
        {
            return new HybridVect(-left / right.value, right.index, right.genealogy);
        }
        return new HybridVect(left / right.value, right.index, right.genealogy);
    }

    public static bool operator ==(HybridVect left, HybridVect right)
    {
        return (left.value == right.value && left.index == right.index && left.genealogy == right.genealogy);
    }
    public static bool operator !=(HybridVect left, HybridVect right)
    {
        return !(left == right);
    }

    public static int DefaultSigfigs
    {
        get
        {
            return Vect.DefaultSigfigs;
        }
        set
        {
            Vect.DefaultSigfigs = value;
        }
    }
    // Convert this Vector to a String.
    public override string ToString()
    {
        string s = this.value.RoundString(DefaultSigfigs);
        if (s.Contains("E"))
        {
            s = "(" + s + ")";
        }
        return s + "e" + ToSubscript(this.index);
    }
    // Convert this Vector to a String, with the provided rules.
    // PrintFlags.MatrixForm and PrintFlags.Reorder are currently not supported,
    // But PrintFlags.ForceSign, PrintFlags.HideZeroes, and PrintFlags.SpaceSign are
    public string ToString(PrintFlags flags)
    {
        string res = "";
        if ((flags & PrintFlags.MatrixForm) != 0)
        {
            throw new NotImplementedException();
        }
        if ((flags & PrintFlags.Reorder) != 0)
        {
            throw new NotImplementedException();
        }
        if ((flags & PrintFlags.HideZeroes) != 0 && this.value == 0)
        {
            return "";
        }
        bool isNegative = this.value < 0;
        string s = (isNegative ? -this.value : this.value).RoundString(Vect.DefaultSigfigs);
        if (s.Contains("E"))
        {
            s = "(" + s + ")";
        }
        if ((flags & PrintFlags.ForceSign) != 0 && (this.value >= 0 || s.Contains("E")))
        {
            res += "+";
        }
        else if (this.value < 0)
        {
            res += "-";
        }

        if ((flags & PrintFlags.SpaceSign) != 0)
        {
            res = " " + res + " ";
        }
        res += s + "e" + ToSubscript(this.index);
        return res;
    }

    // Convert a uint to a subscript string.
    private static string ToSubscript(uint num)
    {
        if (num == 0)
        {
            return "₀";
        }
        string res = "";
        while (num > 0)
        {
            char val = (char)(num % 10);
            num /= 10;
            res = ((char)(val + '₀')).ToString() + res;
        }
        return res;
    }
}

public struct HybridTessarine
{
    // Create a Hybrid Tessarine from a list of Doubles
    public HybridTessarine(uint genealogy, params double[] vects)
    {
        HybridVect[] arr = new HybridVect[vects.Length.CeilDimension()];
        for (int i = 0; i < vects.Length; i++)
        {
            arr[i] = new HybridVect(vects[i], i, genealogy);
        }
        for (int i = vects.Length; i < arr.Length; i++)
        {
            arr[i] = new HybridVect(0, i, genealogy);
        }
        this.vects = arr;
        this.dim = (uint)arr.Length;
        this._genealogy = genealogy;
    }
    // Create a Tessarine from a list of Vectors
    public HybridTessarine(params HybridVect[] vects)
    {
        if (!vects.IsAllSameGenealogy()) throw new InvalidOperationException("Cannot create a Hybrid Tessarine from vectors of differing genealogies");
        this._genealogy = vects[0].genealogy;
        this.vects = null;
        this.dim = 2;

        HybridVect[] arr = CreateVectArray(vects.GetHighestDimension(), this.genealogy);
        for (int i = 0; i < vects.Length; i++)
        {
            HybridVect v = vects[i];
            arr[v.index].value += v.value;
        }
        this.vectors = arr;
    }

    public HybridTessarine(uint genealogy, Tessarine proper)
    {
        HybridVect[] v = new HybridVect[proper.dimensions];
        for (int i = 0; i < v.Length; i++)
        {
            v[i] = new HybridVect(proper[i], genealogy);
        }
        this.vects = CreateVectArray(v);
        this.dim = proper.dimensions;
        this._genealogy = genealogy;
    }
    public bool IsZeroDivisor()
    {
        return ((Tessarine)this).IsZeroDivisor();
    }
    public static HybridVect FirstRotator(uint genealogy)
    {
        return new HybridVect(1, (1 << genealogy.LowestSetBitIndex()), genealogy);
    }

    public double Abs()
    {
        if (this.dimensions == 2)
        {
            double a = this[0].value;
            double b = this[1].value;
            if ((this.genealogy & 1) == 0)
            {
                return Math.Sqrt(Math.Abs(a * a - b * b));
            }
            else
            {
                return Math.Sqrt(a * a + b * b);
            }
        }

        if ((this.genealogy & 1) == 0 && this.genealogy != 0)
        {
            return new HybridTessarine(this.genealogy.SwapBits(this.genealogy.LowestSetBitIndex(), 0), (Tessarine)this).Abs();
        }

        HybridTessarine red = this.Reduce();

        HybridTessarine L = red.left;
        HybridTessarine R = red.right;

        if ((this.genealogy & (this.dimensions >> 1)) != 0)
        {
            R *= FirstRotator(this.genealogy);
        }

        return Math.Sqrt((L + R).Abs() * (L - R).Abs());
    }

    // Get the sum of absolute-values of each vector.
    // Useful for finding the difference between two Tessarines,
    // Particularly for finding the computational error in certain operations
    public double AbsMagnitude()
    {
        double d = 0;
        for (int i = 0; i < this.vects.Length; i++)
        {
            d += Math.Abs(this.vects[i].value);
        }
        return d;
    }

    // Convert this Hybrid Tessarine to a String (Plain text)
    // I have yet to implement the inverse, but I think it's better not to have it
    public override string ToString()
    {
        string res = this[0].ToString();
        for (int i = 1; i < this.dimensions; i++)
        {
            res += this[i].ToString(PrintFlags.ForceSign | PrintFlags.SpaceSign);
        }
        return "[" + Convert.ToString(this.genealogy & (this.dimensions - 1), 2) + "] " + res;
    }
    // Returns 'true' if every term in this Tessarine is zero; Otherwise, 'false'
    public bool IsZero()
    {
        for (int i = 0; i < this.vects.Length; i++)
        {
            if (this.vects[i].value != 0)
            {
                return false;
            }
        }
        return true;
    }

    // Creates a zeroed Vect array of the given length (Rounded *up* to the closest suitable dimension)
    public static HybridVect[] CreateVectArray(int length, uint genealogy)
    {
        return CreateVectArray((uint)length, genealogy);
    }
    // Same as above but better
    public static HybridVect[] CreateVectArray(uint length, uint genealogy)
    {
        HybridVect[] arr = new HybridVect[length.CeilDimension()];
        for (int i = 0; i < arr.Length; i++)
        {
            arr[i] = new HybridVect(0, i, genealogy);
        }
        return arr;
    }
    // Takes in a Vect array, fixes the ordering (And adds same-index Vectors together),
    // And returns a new Vect[] array with the smallest suitable dimensions
    public static HybridVect[] CreateVectArray(HybridVect[] arr)
    {
        if (!arr.IsAllSameGenealogy()) throw new InvalidOperationException("Cannot create a vector array from vectors of differing genealogies");
        HybridVect[] res = CreateVectArray(arr.GetHighestDimension(), arr[0].genealogy);
        for (int i = 0; i < arr.Length; i++)
        {
            res[arr[i].index].value += arr[i].value;
        }
        return res;
    }
    // Fixes the ordering of the terms in this Hybrid Tessarine if they are fucked up for some reason (You have to save the result to a variable)
    public HybridTessarine FixOrdering()
    {
        HybridVect[] arr = new HybridVect[this.vectors.GetHighestDimension()];
        for (int i = 0; i < arr.Length; i++)
        {
            arr[i] = new HybridVect(0, i, this.genealogy);
        }
        for (int i = 0; i < this.vects.Length; i++)
        {
            HybridVect v = this.vects[i];
            HybridVect v1 = arr[v.index];
            v1.value += v.value;
            arr[v.index] = v1;
        }
        this.vectors = arr;
        return this;
    }
    // Creates a copy of this Hybrid Tessarine
    public HybridTessarine Copy()
    {
        return new HybridTessarine(HybridTessarine.CreateVectArray(this.vectors));
    }
    // Reduces this Hybrid Tessarine to the smallest-dimensional Hybrid Tessarine that represents the exact same value.
    public HybridTessarine Reduce()
    {
        if (this.dimensions == 2) return this.Copy();
        HybridVect[] v = this.vectors;
        uint i = this.dimensions >> 1;
        for (; i < v.Length; i++)
        {
            if (v[i].value != 0)
            {
                return this.Copy();
            }
        }
        return this.left.Reduce();
    }
    // Resizes this Hybrid Tessarine to the dimensions given (Ceiled/Rounded up to closest suitable dimension count)
    public HybridTessarine ToDimension(uint dimension)
    {
        dimension = dimension.CeilDimension();
        HybridVect[] newVects = CreateVectArray(dimension, this.genealogy);
        uint min = dimension;
        if (this.dimensions < min) min = this.dimensions;

        for (int i = 0; i < min; i++)
        {
            newVects[i] = this[i];
        }
        this.vectors = newVects;
        return this;
    }

    public static HybridTessarine Random(uint genealogy, uint dimensions)
    {
        return Random(genealogy, dimensions, 0, 1);
    }
    public static HybridTessarine Random(uint genealogy, int dimensions)
    {
        return Random(genealogy, (uint)dimensions, 0, 1);
    }
    public static HybridTessarine Random(uint genealogy, int dimensions, double min, double max)
    {
        return Random(genealogy, (uint)dimensions, min, max);
    }
    public static HybridTessarine Random(uint genealogy, uint dimensions, double min, double max)
    {
        HybridVect[] v = new HybridVect[dimensions];
        double range = max - min;
        for (int i = 0; i < dimensions; i++)
        {
            v[i] = new HybridVect((Global.rand.NextDouble() * range) + min, i, genealogy);
        }
        return new HybridTessarine(v);
    }


    // ##############


    // Multiply two Tessarines
    public static HybridTessarine operator *(HybridTessarine left, HybridTessarine right)
    {
        if (left.genealogy != right.genealogy) throw new InvalidOperationException("Cannot multiply Tessarines of differing genealogies");
        HybridVect[] mulBV = right.vectors;
        HybridVect[] addVals = new HybridVect[mulBV.Length];
        HybridTessarine t = new HybridTessarine(left.genealogy);
        for (int i = 0; i < left.dimensions; i++)
        {
            if (left[i].value == 0) continue;
            for (int j = 0; j < mulBV.Length; j++)
            {
                addVals[j] = left[i] * mulBV[j];
            }
            t += addVals;
        }
        return t;
    }

    // Multiply a Tessarine by a Vector
    public static HybridTessarine operator *(HybridTessarine left, HybridVect right)
    {
        if (left.genealogy != right.genealogy) throw new InvalidOperationException("Cannot multiply differing genealogies");
        if (left.dimensions <= right.index)
        {
            left = left.ToDimension(Math.Max(left.dimensions, (right.index + 1).CeilDimension()));
        }
        HybridVect[] vectArr = left.vectors;
        HybridVect[] resArr = CreateVectArray(vectArr.Length, left.genealogy);
        HybridTessarine t = new HybridTessarine(left.genealogy);
        for (int i = 0; i < vectArr.Length; i++)
        {
            HybridVect v = vectArr[i] * right;
            resArr[v.index] = v;
        }
        t.vectors = resArr;
        return t;
    }

    // Multiply a Vector by a Tessarine
    public static HybridTessarine operator *(HybridVect left, HybridTessarine right)
    {
        return right * left;
    }

    // Multiply a Tessarine by an array of Vectors
    public static HybridTessarine operator *(HybridTessarine left, HybridVect[] right)
    {
        if (!right.IsAllSameGenealogy() || left.genealogy != right[0].genealogy) throw new InvalidOperationException("Cannot multiply Tessarines of differing genealogies");
        HybridTessarine res = new HybridTessarine(left.genealogy);
        for (int i = 0; i < right.Length; i++)
        {
            res += left * right[i];
        }
        return res;
    }

    // Multiply a Vector array by a Tessarine
    public static HybridTessarine operator *(HybridVect[] left, HybridTessarine right)
    {
        return right * left;
    }

    // Add a Tessarine and a double
    public static HybridTessarine operator +(HybridTessarine left, double right)
    {
        HybridTessarine t = left.Copy();
        HybridVect v = t[0];
        v.value += right;
        t[0] = v;
        return t;
    }

    // Same as above
    public static HybridTessarine operator +(double left, HybridTessarine right)
    {
        return right + left;
    }

    // Subtract a double from a Tessarine
    public static HybridTessarine operator -(HybridTessarine left, double right)
    {
        HybridTessarine t = left.Copy();
        HybridVect v = t[0];
        v.value -= right;
        t[0] = v;
        return t;
    }

    // Subtract a Tessarine from a double
    public static HybridTessarine operator -(double left, HybridTessarine right)
    {
        HybridTessarine t = -right.Copy();
        HybridVect v = t[0];
        v.value = left + v.value;
        t[0] = v;
        return t;
    }

    // Multiply a Tessarine by a double
    public static HybridTessarine operator *(HybridTessarine left, double right)
    {
        HybridVect[] v = left.vectors.Copy();
        for (int i = 0; i < v.Length; i++)
        {
            v[i] *= right;
        }
        return new HybridTessarine(v);
    }

    // Same as above
    public static HybridTessarine operator *(double left, HybridTessarine right)
    {
        return right * left;
    }

    // Divide a Tessarine by a double;
    public static HybridTessarine operator /(HybridTessarine left, double right)
    {
        HybridVect[] v = left.vectors.Copy();
        right = 1 / right;
        for (int i = 0; i < v.Length; i++)
        {
            v[i] *= right;
        }
        return new HybridTessarine(v);
    }

    // Divide a double by a Tessarine
    public static HybridTessarine operator /(double left, HybridTessarine right)
    {
        if (left == 1)
        {
            return TessMath.exp(-TessMath.ln(right));
        }
        else if (left == -1)
        {
            return -TessMath.exp(-TessMath.ln(right));
        }
        return left * TessMath.exp(-TessMath.ln(right));
    }

    // Divide a Vector by a Tessarine
    public static HybridTessarine operator /(HybridVect left, HybridTessarine right)
    {
        if (left.genealogy != right.genealogy) throw new InvalidOperationException("Cannot create a Hybrid Tessarine from vectors of differing genealogies");
        return left * TessMath.exp(-TessMath.ln(right));
    }

    // Divide a Tessarine by a Vector
    public static HybridTessarine operator /(HybridTessarine left, HybridVect right)
    {
        if (left.genealogy != right.genealogy) throw new InvalidOperationException("Cannot create a Hybrid Tessarine from vectors of differing genealogies");
        HybridVect[] v = left.vectors;
        for (int i = 0; i < v.Length; i++)
        {
            v[i] = v[i] / right;
        }
        return new HybridTessarine(v);
    }

    // Divide a Tessarine by a Tessarine
    public static HybridTessarine operator /(HybridTessarine left, HybridTessarine right)
    {
        if (left.genealogy != right.genealogy) throw new InvalidOperationException("Cannot create a Hybrid Tessarine from vectors of differing genealogies");
        return TessMath.div(left, right);
    }


    // Add a Tessarine and a Vector array
    public static HybridTessarine operator +(HybridTessarine left, HybridVect[] right)
    {
        if (!right.IsAllSameGenealogy() || left.genealogy != right[0].genealogy) throw new InvalidOperationException("Cannot multiply Tessarines of differing genealogies");
        HybridTessarine t = left.Copy();
        t = t.ToDimension(Math.Max(left.dimensions, right.GetHighestDimension()));

        for (int i = 0; i < right.Length; i++)
        {
            HybridVect r = right[i];
            HybridVect v = t[r.index];
            v.value += r.value;
            t[r.index] = v;
        }
        return t;
    }

    // Same as above
    public static HybridTessarine operator +(HybridVect[] left, HybridTessarine right)
    {
        return right + left;
    }

    // Add a Tessarine to another Tessarine
    public static HybridTessarine operator +(HybridTessarine left, HybridTessarine right)
    {
        if (left.genealogy != right.genealogy)
        {
            throw new InvalidOperationException("Cannot create a Hybrid Tessarine from vectors of differing genealogies");
        }
        return left + right.vectors;
    }

    // Just copies the Tessarine value
    public static HybridTessarine operator +(HybridTessarine val)
    {
        return val.Copy();
    }

    // Subtract a Vector array from a Tessarine
    public static HybridTessarine operator -(HybridTessarine left, HybridVect[] right)
    {
        if (!right.IsAllSameGenealogy() || left.genealogy != right[0].genealogy) throw new InvalidOperationException("Cannot multiply Tessarines of differing genealogies");
        HybridTessarine t = left.Copy();
        t = t.ToDimension(Math.Max(left.dimensions, right.GetHighestDimension()));

        for (int i = 0; i < right.Length; i++)
        {
            HybridVect r = right[i];
            HybridVect v = t[r.index];
            v.value -= r.value;
            t[r.index] = v;
        }
        return t;
    }

    // Subtract a Tessarine from a Vector array
    public static HybridTessarine operator -(HybridVect[] left, HybridTessarine right)
    {
        if (!left.IsAllSameGenealogy() || right.genealogy != left[0].genealogy) throw new InvalidOperationException("Cannot multiply Tessarines of differing genealogies");
        HybridTessarine t = right.Copy();
        t = t.ToDimension(Math.Max(right.dimensions, left.GetHighestDimension()));

        for (int i = 0; i < left.Length; i++)
        {
            HybridVect r = left[i];
            HybridVect v = t[r.index];
            v.value = r.value - v.value;
            t[r.index] = v;
        }
        return t;
    }

    // Subtract a Tessarine from a Tessarine
    public static HybridTessarine operator -(HybridTessarine left, HybridTessarine right)
    {
        if (left.genealogy != right.genealogy) throw new InvalidOperationException("Cannot create a Hybrid Tessarine from vectors of differing genealogies");
        return left - right.vectors;
    }

    // Negate this Tessarine
    public static HybridTessarine operator -(HybridTessarine val)
    {
        HybridVect[] vects = val.vectors.Copy();
        for (int i = 0; i < vects.Length; i++)
        {
            vects[i] = -vects[i];
        }
        return new HybridTessarine(vects);
    }
    // ##############

    // Exponentiate a Tessarine by another Tessarine
    public static HybridTessarine operator ^(HybridTessarine left, HybridTessarine right)
    {
        return TessMath.pow(left, right);
    }

    // Take the sqrt of this Tessarine
    public HybridTessarine sqrt()
    {
        return TessMath.sqrt(this);
    }

    // Square this Tessarine
    public HybridTessarine sq()
    {
        return this * this;
    }

    // Get the reciprocal of this Tessarine
    public HybridTessarine reciprocal()
    {
        return 1 / this;
    }

    // Return the exp (e^x) of this Tessarine
    public HybridTessarine exp()
    {
        return TessMath.exp(this);
    }

    // Compute the natural log of this Hybrid Tessarine
    public HybridTessarine ln()
    {
        return TessMath.ln(this);
    }

    // Compute the conjugate of this Hybrid Tessarine
    public Tessarine conj()
    {
        return TessMath.conj(this);
    }

    public static implicit operator HybridVect[](HybridTessarine val)
    {
        return val.vectors.Copy();
    }
    public static implicit operator double[](HybridTessarine val)
    {
        double[] res = new double[val.dimensions];
        for (int i = 0; i < res.Length; i++)
        {
            res[i] = val[i].value;
        }
        return res;
    }
    public static implicit operator Tessarine(HybridTessarine val)
    {
        Vect[] v = new Vect[val.dimensions];
        for (int i = 0; i < v.Length; i++)
        {
            v[i] = val[i];
        }
        return new Tessarine(v);
    }
    public static implicit operator HybridTessarine(HybridVect v)
    {
        return new HybridTessarine(v);
    }
    public static implicit operator HybridTessarine(HybridVect[] v)
    {
        return new HybridTessarine(v);
    }
    public static bool operator true(HybridTessarine t)
    {
        return !t.IsZero();
    }

    // Is this Tessarine zero?
    public static bool operator false(HybridTessarine t)
    {
        return t.IsZero();
    }
    // Get or set the vector at the given index
    // This is where you can fuck up the ordering on the 'set' part
    public HybridVect this[int index]
    {
        get
        {
            if (index >= this.vectors.Length)
            {
                return new HybridVect(0, index, this.genealogy);
            }
            return this.vectors[index];
        }
        set
        {
            if (index >= this.vectors.Length)
            {
                HybridVect[] arr = CreateVectArray(index, this.genealogy);
                for (int i = 0; i < this.vectors.Length; i++)
                {
                    HybridVect v = arr[this.vectors[i].index];
                    v.value += this.vectors[i].value;
                    arr[this.vectors[i].index] = v;
                }
                this.vectors = arr;
            }
            this.vectors[index] = value;
        }
    }

    // Same as above but no negative indexes
    public HybridVect this[uint index]
    {
        get
        {
            if (index >= this.vectors.Length)
            {
                return new HybridVect(0, index, this.genealogy);
            }
            return this.vectors[(int)index];
        }
        set
        {
            if (index >= this.vectors.Length)
            {
                HybridVect[] arr = CreateVectArray(index, this.genealogy);
                for (int i = 0; i < this.vectors.Length; i++)
                {
                    HybridVect v = arr[this.vectors[i].index];
                    v.value += this.vectors[i].value;
                    arr[this.vectors[i].index] = v;
                }
                this.vectors = arr;
            }
            this.vectors[(int)index] = value;
        }
    }

    // Get the left half of this Tessarine
    public HybridTessarine left
    {
        get
        {
            HybridVect[] v = new HybridVect[this.dimensions >> 1];
            for (int i = 0; i < v.Length; i++)
            {
                v[i] = this.vectors[i];
            }
            return new HybridTessarine(v);
        }
    }

    // Get the right half of this Tessarine
    public HybridTessarine right
    {
        get
        {
            uint p = this.dimensions >> 1;
            HybridVect[] v = new HybridVect[p];
            for (int i = 0; i < v.Length; i++)
            {
                HybridVect v1 = this.vectors[v.Length + i];
                v1.index ^= p;
                v[i] = v1;
            }
            return new HybridTessarine(v);
        }
    }
    // Get or set the vectors in this Hybrid Tessarine
    public HybridVect[] vectors
    {
        get
        {
            return this.vects;
        }
        set
        {
            this.vects = CreateVectArray(value);
            this.dimensions = (uint)this.vects.Length;
        }
    }

    // Get the dimensions of this Hybrid Tessarine
    public uint dimensions
    {
        get
        {
            return dim;
        }
        private set
        {
            dim = value;
        }
    }
    public uint genealogy
    {
        get
        {
            return this._genealogy;
        }
        set
        {
            this.vects = new HybridTessarine(value, (Tessarine)this).vects;
            this._genealogy = value;
        }
    }
    private uint dim;
    private HybridVect[] vects;
    private uint _genealogy;

}
