using System;
using System.Collections.Generic;
using System.Collections;
using System.Linq;
using System.Text;
using System.Numerics;
using System.Security.Cryptography;
using Emil.GMP;


namespace Pollard
{
    class Program
    {
        static int rowIndex = -1;
        static BigInteger GCD;
        static bool found;
        static bool copyMatrix;
        static int indexCopy = 0;
        static int rCopy = 0;
        static int uCopy = 0;
        static int BabyN = 3;// Set the BabyN for babyMethod() even is 2n + 2 odd is 2n + 1
        static bool iNumIsJacobiPrime = true;
        static string[] lines = System.IO.File.ReadAllLines(@"C:\Cryptography\Pollard\Pollard\primes.txt");    
        static BigInt[] xValues;
        static BigInt[] result;
        static int size =  lines.Count()-2;
        static BigInt[] BSmoothCopy = new BigInt[size];
        static UInt64[] tempArray;
        static UInt64[] JacobiPrime;
        static BigInt[,] m = new BigInt[10000, 10000];
        static BigInt[,] BSmoothMatrix = new BigInt[10000, 10000];
        static BigInt[] BSmooth = new BigInt[size];
        static int BSmoothLoop = 0;
        private static double gcd(double a, Int64 b)
        {
            Int64 t;

            // Ensure B > A
            if (Convert.ToInt64(a) > b)
            {
                t = b;
                Convert.ToDouble(b);
                b = Convert.ToInt32(a);
                a = Convert.ToDouble(t);
            }

            // Find 
            while (b != 0)
            {
                t = Convert.ToInt64(a) % b;
                a = (double)b;
                b = t;
            }

            return a;
        }

        static BigInteger SimpleFactor()
        {
            var rng = new RNGCryptoServiceProvider();
            byte[] bytes = new byte[8];
            rng.GetBytes(bytes);
            BigInteger a = new BigInteger(bytes);
            BigInteger b = new BigInteger(bytes);


            Console.WriteLine("Please enter your integer: ");
            a = BigInteger.Parse(Console.ReadLine());
            for (b = 1; b <= a; b++)
            {
                if (a % b == 0)
                {
                    Console.WriteLine(b + " is a factor of " + a);
                }
            }

            return 1;



        }
/*
        static BigInteger Pollard(BigInteger p)
        {
            int COUNT = 2;
            Random rand = new Random();
            // Generate and display 5 random byte (integer) values. 
            byte[] bytes = new byte[4];
            rand.NextBytes(bytes);

            BigInteger[] BaseArray = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113 };
            BigInteger[] L = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
            BigInteger[] expBaseArray = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
            BigInteger[] checkGCD = { };
            L[0] = (BigInteger)(BigInteger.Log(p) / BigInteger.Log(2));

            for (int loop = 0; loop < 29; loop++)
            {
                L[loop] = (BigInteger)(BigInteger.Log(p) / BigInteger.Log(BaseArray[loop]));
                L[loop] = (BigInteger)Math.Floor((double)L[loop]);

                expBaseArray[loop] = BigInteger.ModPow(BaseArray[loop], L[loop], p);


            }
            BigInteger[] exp = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
            BigInteger[] checkGcd = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
            do
            {
                Console.WriteLine("SEARCHING  "+ COUNT);
                COUNT++;
                int randomA = rand.Next(3, 1000000);
                checkGcd[0] = BigInteger.ModPow(randomA, expBaseArray[0], p);


                exp[0] = BigInteger.ModPow(BaseArray[1], L[1], p);

                try
                {
                    GCD = BigInteger.GreatestCommonDivisor((BigInteger)checkGcd[0] - 1, p);
                }
                catch (ArgumentOutOfRangeException e)
                {
                    Console.WriteLine("Unable to calculate the greatest common divisor:");
                    Console.WriteLine("   {0} is an invalid value for {1}",
                                      e.ActualValue, e.ParamName);

                }


                //////
                if ((GCD) == 1)
                {
                    result[0] = BigInteger.ModPow(checkGcd[0], BigInteger.ModPow(BaseArray[1], L[1], p), p);

                }

                else
                {

                    BigInteger factor = BigInteger.GreatestCommonDivisor((BigInteger)checkGcd[0] - 1, p);
                    Console.WriteLine(factor);
                    found = true;

                }

                //////

                //////
                if (BigInteger.GreatestCommonDivisor((BigInteger)result[0] - 1, p) == 1)
                {
                    exp[1] = BigInteger.ModPow(BaseArray[2], L[2], p);


                    result[1] = BigInteger.ModPow(result[0], (BigInteger)exp[1], p);

                }
                else
                {

                    BigInteger factor = BigInteger.GreatestCommonDivisor((BigInteger)result[0] - 1, p);

                    Console.WriteLine(factor);
                    found = true;
                }
                //////

                /////
                if (BigInteger.GreatestCommonDivisor((BigInteger)result[1] - 1, p) == 1)
                {
                    exp[2] = BigInteger.ModPow(BaseArray[3], L[3], p);


                    result[2] = BigInteger.ModPow(result[1], (BigInteger)exp[2], p);

                }

                else
                {


                    BigInteger factor = BigInteger.GreatestCommonDivisor((BigInteger)result[1] - 1, p);
                    Console.WriteLine(factor);
                    found = true;
                }
                //////



                /////
                if (BigInteger.GreatestCommonDivisor((BigInteger)result[2] - 1, p) == 1)
                {
                    exp[3] = BigInteger.ModPow(BaseArray[4], L[4], p);


                    result[3] = BigInteger.ModPow(result[2], (BigInteger)exp[3], p);

                }

                else
                {


                    BigInteger factor = BigInteger.GreatestCommonDivisor((BigInteger)result[2] - 1, p);
                    Console.WriteLine(factor);
                    found = true;
                }
                //////


                /////
                if (BigInteger.GreatestCommonDivisor((BigInteger)result[3] - 1, p) == 1)
                {
                    exp[4] = BigInteger.ModPow(BaseArray[5], L[5], p);


                    result[4] = BigInteger.ModPow(result[3], (BigInteger)exp[4], p);

                }

                else
                {


                    BigInteger factor = BigInteger.GreatestCommonDivisor((BigInteger)result[3] - 1, p);
                    Console.WriteLine(factor);
                    found = true;
                }
                //////


                /////
                if (BigInteger.GreatestCommonDivisor((BigInteger)result[4] - 1, p) == 1)
                {
                    exp[5] = BigInteger.ModPow(BaseArray[6], L[6], p);
                    result[5] = BigInteger.ModPow(result[4], (BigInteger)exp[5], p);

                }

                else
                {


                    BigInteger factor = BigInteger.GreatestCommonDivisor((BigInteger)result[4] - 1, p);
                    Console.WriteLine(factor);

                }
                //////

                /////
                if (BigInteger.GreatestCommonDivisor((BigInteger)result[5] - 1, p) == 1)
                {
                    exp[6] = BigInteger.ModPow(BaseArray[7], L[7], p);


                    result[6] = BigInteger.ModPow(result[5], (BigInteger)exp[6], p);

                }

                else
                {

                    BigInteger factor = BigInteger.GreatestCommonDivisor((BigInteger)result[5] - 1, p);
                    Console.WriteLine(factor);
                    found = true;
                }
                //////
                /////
                if (BigInteger.GreatestCommonDivisor((BigInteger)result[6] - 1, p) == 1)
                {
                    exp[7] = BigInteger.ModPow(BaseArray[8], L[8], p);
                    result[6] = BigInteger.ModPow(result[6], (BigInteger)exp[7], p);

                }

                else
                {


                    BigInteger factor = BigInteger.GreatestCommonDivisor((BigInteger)result[6] - 1, p);
                    Console.WriteLine(factor);
                    found = true;
                }


                for (int x = 1; x < 16; x++)
                {

                    if (BigInteger.GreatestCommonDivisor((BigInteger)result[6 + x] - 1, p) == 1)
                    {
                        exp[7 + x] = BigInteger.ModPow(BaseArray[8 + x], L[8 + x], p);


                        result[6 + x] = BigInteger.ModPow(result[6 + x], (BigInteger)exp[7 + x], p);


                    }

                    else
                    {

                        BigInteger factor = BigInteger.GreatestCommonDivisor((BigInteger)result[6 + 1] - 1, p);
                        Console.WriteLine(factor);
                        found = true;
                    }
                }

            }
            while (found == false);
            return 1;

        }
 */

        static BigInt primeFactors(BigInt N, BigInt p)
        {
           // BigInt nSqrt = BabylonianMethod(N);

            p = new BigInt(p);
            BigInt[] tempArray;
           
            
            tempArray = new BigInt[size];
        
         
            int countExp = 0;
            int exponent = 1;
            bool exit = false;
            int columnIndex = 0;
            rowIndex++;
            BigInt orignalN = N;
       



            for (BigInt i = 2; i <= N; i += 1)
            {
                
                if (N % i == 0)
                {
                    while (N % i == 0)
                    {
                        exit = false;
                        countExp = 0;
                        columnIndex = 0;

                        N /= i;
                        if (i.IsProbablyPrimeRabinMiller(200))
                        { 
                            if((rowIndex + p+1)==336)
                            {
                                Console.Write(" ");

                        }
                            Console.WriteLine("Factors of " + orignalN + "   " + i + "  RowIndex " + (rowIndex + p+1));
                            Console.WriteLine(" ");
                            //Console.WriteLine(JacobiPrime.Count());

                        }
                        while (exit == false)
                        {
                            //if i matches JacobiPrime[countExp]
                           
                            if (i == JacobiPrime[columnIndex])
                            {
                               // Console.WriteLine("Hello " + i + " is in JacobiPrime[" + countExp + "]");
                                  m[rowIndex, columnIndex] += exponent;
                                //Console.WriteLine(m[rowIndex, columnIndex]);
                                exit = true;
                               
                            }
                            //else store a 0 in the position increment  countExp++
                            else
                            {
                               
                                if (countExp < 150)
                                {
                                    columnIndex++;
                                    countExp++;
                                   
                                }
                                else {
                                    if (N % i == 0 )
                                    {
                                        iNumIsJacobiPrime = false;
                                    }
                                    
                                    
                                    exit = true;  }
                            }
                        }
                        
                        if (i >= N && i < 200 &&  iNumIsJacobiPrime == true)//if it is BSmooth
                        {
                            BSmooth[BSmoothLoop] = orignalN; //Store BSmooth Numbers
                            BSmoothLoop++;
                         
                          
                        }
                    


                     }
                    columnIndex++;
                }

              
            }

            return 1;
          
           
        }

        static BigInt BabylonianMethod(BigInt SquareNumber)
        {
             
            BigInt result = new BigInt(SquareNumber);
 
            BigInt[] values;
            values = new BigInt[1000000];
           
            int loop = 1;
            values[loop - 1] = 2* BigInt.Power(10, BabyN);
            bool exit = false;
            do {
                values[loop] = (values[loop - 1] + SquareNumber / values[loop - 1]) / 2;
                //Console.WriteLine(values[loop]);
                if (values[loop] == values[loop - 1])
                { 
                    exit = true; }
                loop++;
            } while (!exit);

            result = (values[loop - 1]);

            return result;

        }

        static int  FindFactorBase(BigInt RSAModulus)
        {
            tempArray = new UInt64[size];
             
            JacobiPrime = new UInt64[size];
            for (int LineLoop = 0; LineLoop < lines.Count()-2; LineLoop++)
            {

               tempArray[LineLoop] =  Convert.ToUInt64(lines[LineLoop]);
               //Console.WriteLine(lines[LineLoop]);
            }
            UInt64 jacobiLoop = 1;
            JacobiPrime[0] = 2;
            Console.WriteLine("Jacobi Prime " + JacobiPrime[0]);
            for (int primeLoop = 1; primeLoop <151; primeLoop++)
            {
                if (BigInt.JacobiSymbol(RSAModulus,tempArray[primeLoop])==1)
                {
                    JacobiPrime[jacobiLoop] = Convert.ToUInt64(tempArray[primeLoop]);
                 Console.WriteLine("Jacobi Prime " + JacobiPrime[jacobiLoop]);
                    jacobiLoop++;
                }
            }
            return 1;
        }

        static BigInt FindxValues(BigInt aproxSqRoot, BigInt RSAmodulus)
        {
          
            Console.WriteLine("Finding x Values.... ");
            xValues = new BigInt[100000];
            
            for(int xValuesLoop =0; xValuesLoop <250; xValuesLoop++){
             
            xValues[xValuesLoop] = aproxSqRoot + xValuesLoop;
            xValues[xValuesLoop] =( xValues[xValuesLoop] * xValues[xValuesLoop]) % RSAmodulus;
             
        //    Console.WriteLine(xValues[xValuesLoop]);
            }


            return 1;

        }
         

        static void Main(string[] args)
        { 
             //BigInteger factor     = Pollard(1288005205276048901);
             //BigInt numberToFactor = new BigInt("1522605027922533360535618378132637429718068114961380688657908494580122963258952897654000350692006139");
             //BigInt primeFactors   = new BigInt("43423");
             //BigInt factorQS       = QS(primeFactors);
             //Console.WriteLine(numberToFactor);
            for (int mLoopCol = 0; mLoopCol < 10000; mLoopCol++)
            {
                for (int mLoopRow = 0; mLoopRow < 10000; mLoopRow++)
                {
                    m[mLoopRow, mLoopCol] = 0;


                }


            }
            BigInt numberToFactor = new BigInt("43423"); // Our RSA modulus to Factor into two Prime Numbers
            BigInt print = BabylonianMethod(numberToFactor);  // Estimate the Sqaure Root of Large Compostite Numbers
            Console.WriteLine(print);
            Console.WriteLine(numberToFactor.Sqrt());
            FindFactorBase(numberToFactor);
            FindxValues(print, numberToFactor);
            for (int primeFactorsLoop = 1; primeFactorsLoop < 250; primeFactorsLoop++)
            {
                primeFactors(xValues[primeFactorsLoop],print);
            }
           
            Console.WriteLine("-----------------------");
            Console.WriteLine("-----------------------");
            Console.WriteLine("                       ");
            
            for(int r = 0;r<150;r++){
            //Console.Write("\n");
            for (int y = 0; y < size; y++)
            {
                 
                if (xValues[r + 1] == BSmooth[y])
                {
                   // Console.Write(xValues[r + 1] + " ");
                    BSmoothCopy[indexCopy] = xValues[r + 1];
                    indexCopy++;
                   
                    for (int u = 0; u < 150; u++)
                    {
                        BSmoothMatrix[rCopy, u] = m[r, u];
                    }
                    rCopy++;
                }
            }
            
            }
          
            for (int row = 0; row <150; row++)
            {
                Console.Write(BSmoothCopy[row] + " ");
                for (int col = 0; col < 150; col++)
                {
                    Console.Write(BSmoothMatrix[row, col]  );
                }

                Console.Write("\n");
            }

            Console.WriteLine("GCD( N, X+Y )= ");
            Console.Write(BigInteger.GreatestCommonDivisor(43423, (407844343998 + 158037298298880)));

            Console.Read();



        }


    }
}




