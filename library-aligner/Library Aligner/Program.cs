using System;

using System.IO;
using System.Linq;
using System.Collections.Generic;
using System.Xml.Serialization;
using System.Xml;
using System.Collections.Concurrent;
using System.Text;
using System.Threading;
using System.Threading.Tasks;

using Library_Analysis;

namespace Library_Aligner
{
    class Program
    {
        static void Main(string[] args)
        {
            //Accessories.CombineFiles(@"C:\Temp\Virome\");
            try
            {
                string XMLPath;
                if (args.Length == 0) XMLPath = ""; else XMLPath = args[0];
                Library_Aligner.Start_Library_Aligner(XMLPath, true);  //gRNA Library Analysis

                //Below are other (Accessories) programs you can run from here
                bool skip = true;
                if (!skip)
                {
                    Accessories.Only_Join();

                    Accessories.Start_Demux();
                    Accessories.Demux_gRNA_Sequences2Files(); //2018.08.14
                    Accessories.CountUMIsFolder(@"E:\Temp\NGS\NGS\Demux\");

                    Accessories.Start_Protein_Aligner();  //Protein Dynamic Alignment
                    gRNA_OffTarget.Find_OffTargets.A_StartHere();

                    Accessories.Start_STARS_Analysis();   //Depletion Analysis for CRISPR Pools

                    DistanceMatrix.aaaTry();
                }
            }
            catch (Exception E)
            {
                Console.WriteLine("");
                Console.WriteLine(E.Message);
                if (E.InnerException != null) Console.WriteLine(E.InnerException.Message);
                Console.ReadKey();
            }
        }
    }
}
