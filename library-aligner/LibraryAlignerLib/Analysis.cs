using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;

namespace Library_Analysis
{
    public class NGS_Results
    {
        public string OrigPath { get; internal set; }
        public List<SampleClass> Samples { get; internal set; }
        public List<string> Sample_Names { get; internal set; }
        public Dictionary<string, GeneClass> Genes { get; internal set; }
        public int GuideCount { get; internal set; }

        public bool RemovedEmptyGuides { get; internal set; }
        public bool LogNormalized { get; internal set; }

        public NGS_Results(string ResRef_FullPath)
        {
            Init(ResRef_FullPath, "0");
        }

        public NGS_Results(string ResRef_FullPath, string ReferenceIndex)
        {
            Init(ResRef_FullPath, ReferenceIndex);
        }

        private static string[] splitTab = new string[1] { "\t" };

        private void Init(string ResRef_FullPath, string ReferenceIndex)
        {
            GuideCount = 0;
            Genes = new Dictionary<string, GeneClass>();
            RemovedEmptyGuides = false;
            LogNormalized = false;

            OrigPath = ResRef_FullPath;
            StreamReader SR = new StreamReader(ResRef_FullPath);
            List<string> Header = SR.ReadLine().Split(splitTab, StringSplitOptions.RemoveEmptyEntries).ToList();
            List<string> LineArray;
            int Col_RefIndex = 0;
            int Col_Name = 1;
            int Col_StartSamples = 7;
            int Col_CountSamples = Header.Count - Col_StartSamples;
            Init_Samples(Header.GetRange(Col_StartSamples, Col_CountSamples));
            while (!SR.EndOfStream)
            {
                LineArray = SR.ReadLine().Split(splitTab, StringSplitOptions.None).ToList();
                if (LineArray[Col_RefIndex] == ReferenceIndex)
                {
                    Add_ResRefLine(LineArray[Col_Name], LineArray.GetRange(Col_StartSamples, Col_CountSamples));
                    GuideCount++;
                }
            }
            SR.Close();
        }

        public void ExportGuides(string Path)
        {
            StreamWriter SW = new StreamWriter(Path);
            bool first = true;
            foreach (GeneClass Gene in Genes.Values)
            {
                foreach (gRNAClass guide in Gene.Guides.Values)
                {
                    if (first)
                    {
                        SW.WriteLine("Long" + "\t" + "Gene" + "\t" + guide.ExportLine(true, Sample_Names));
                        first = false;
                    }
                    SW.WriteLine("" +"\t" + Gene.Name + "\t" + guide.ExportLine(false, null));
                }
            }
            SW.Close();
        }

        public void RemoveGuidesWithFewReads()
        {
            RemoveGuidesWithFewReads(1);
        }

        public void RemoveGuidesWithFewReads(double MinimumReads)
        {
            if (!RemovedEmptyGuides)
            {
                GuideCount = 0;
                foreach (GeneClass Gene in Genes.Values)
                {
                    Gene.RemoveGuideWithFewerThan(MinimumReads);
                    GuideCount += Gene.Guides.Count;
                }
                RemovedEmptyGuides = true;
            }
        }

        public bool AddData_Subtract(string SampleName_Before, string SampleName_After, string DataKeyName)
        {
            int Col1 = -1; int Col2 = -1;
            if (Sample_Names.Contains(SampleName_Before)) Col1 = Sample_Names.IndexOf(SampleName_Before);
            if (Sample_Names.Contains(SampleName_After)) Col2 = Sample_Names.IndexOf(SampleName_After);
            if (Col1 < 0 || Col1 < 0) return false;
            foreach (GeneClass Gene in Genes.Values)
                foreach (gRNAClass guide in Gene.Guides.Values)
                    guide.AddData(DataKeyName, guide.Reads[Col2] - guide.Reads[Col1]);
            return true;
        }

        public bool AddData_Ranks(string DataKeyNameToRank, string RankName_Within, string RankName_Overall)
        {
            SortedList<double, List<gRNAClass>> sortWithin;
            SortedList<double, List<gRNAClass>> sortOverall = new SortedList<double, List<gRNAClass>>();
            double Val;
            foreach (GeneClass Gene in Genes.Values)
            {
                sortWithin = new SortedList<double, List<gRNAClass>>();
                foreach (gRNAClass guide in Gene.Guides.Values)
                {
                    Val = guide.GetData(DataKeyNameToRank);

                    if (!sortOverall.ContainsKey(Val)) sortOverall.Add(Val, new List<gRNAClass>());
                    sortOverall[Val].Add(guide);

                    if (!sortWithin.ContainsKey(Val)) sortWithin.Add(Val, new List<gRNAClass>());
                    sortWithin[Val].Add(guide);
                }
                //Now assign the within ranks
                for (int i = 0; i < sortWithin.Values.Count; i++)
                    foreach (gRNAClass guide in sortWithin.Values[i])
                        guide.AddData(RankName_Within, i + 1);
            }
            //Now assign the overall ranks
            for (int i = 0; i < sortOverall.Values.Count; i++)
                foreach (gRNAClass guide in sortOverall.Values[i])
                    guide.AddData(RankName_Overall, i + 1);
            return true;
        }

        public void Log2Normalize_PerSample()
        {
            if (!LogNormalized)
            {
                foreach (GeneClass Gene in Genes.Values)
                    foreach (gRNAClass gRNA in Gene.Guides.Values)
                        for (int i = 0; i < Samples.Count; i++)
                        {
                            gRNA.Reads[i] = Math.Log(1 + (1000000 * gRNA.Reads[i] / Samples[i].TotalReads));
                        }

                LogNormalized = true;
            }
        }

        private void Init_Samples(List<string> SampleNames)
        {
            Sample_Names = SampleNames;
            Samples = new List<SampleClass>(SampleNames.Count);
            for (int i = 0; i < SampleNames.Count; i++)
            {
                Samples.Add(new SampleClass(SampleNames[i]));
            }
        }

        private void Add_ResRefLine(string Name, List<string> SampleReads)
        {
            GeneClass tGene;
            gRNAClass tGuide;
            bool Success =
                GetAddGene(Name, out tGene, out tGuide);
            double d;
            if (Success)
            {
                for (int i = 0; i < SampleReads.Count; i++)
                {
                    if (!double.TryParse(SampleReads[i], out d)) d = 0;
                    tGuide.AddReads(d);
                    Samples[i].AddReads(d);
                }
            }
        }

        private bool GetAddGene(string Name, out GeneClass Gene, out gRNAClass gRNA)
        {
            Gene = null;
            gRNA = null;
            string tGene, tgRNA, tgRNA_Mod;
            bool Success = ParseName(Name, out tGene, out tgRNA, out tgRNA_Mod);
            if (Success)
            {
                if (!Genes.ContainsKey(tGene)) Genes.Add(tGene, new GeneClass(tGene));
                Gene = Genes[tGene];
                gRNA = Gene.AddOrGetGuide(tgRNA, tgRNA_Mod);
                return true;
            }
            return false;
        }

        private static bool ParseName(string Name, out string GeneName, out string gRNAName, out string PostgRNA)
        {
            Regex regex = new Regex("(.+)(?: |_| G|_G| GRNA|_GRNA|GRNA)(\\d+)(.*)");
            Match M = regex.Match(Name.Trim().ToUpper());
            if (M.Success)
            {
                if (M.Groups.Count == 4)
                {
                    GeneName = M.Groups[1].Value;
                    gRNAName = M.Groups[2].Value;
                    PostgRNA = M.Groups[3].Value;
                    return true;
                }
            }
            GeneName = "";
            gRNAName = "";
            PostgRNA = "";
            return false;
        }
    }

    public class STARS_Analysis
    {
        public string SampleBefore { get; private set; }
        public string SampleAfter { get; private set; }
        public NGS_Results NGSResults { get; private set; }

        public STARS_Analysis(NGS_Results Results, string SampleBefore, string SampleAfter)
        {
            Init(Results, SampleBefore, SampleAfter);
        }

        public STARS_Analysis(string ResRef_FullPath, string ReferenceIndex, string SampleBefore, string SampleAfter)
        {
            NGS_Results tNGSResults = new NGS_Results(ResRef_FullPath, ReferenceIndex);
            Init(tNGSResults, SampleBefore, SampleAfter);
        }

        private void Init(NGS_Results Results, string SampleBefore, string SampleAfter)
        {
            NGSResults = Results;
            this.SampleAfter = SampleAfter;
            this.SampleBefore = SampleBefore;
        }

        private static string nDiff = "LogN_Diff";
        private static string nRankWI = "Rank_Within";
        private static string nRankOA = "Rank_Overall";
        private static string nSTARS = "STARS_Score";
        private static string nSTARSRank = "STARS_Rank";

        public void RunAnalysis()
        {
            NGSResults.RemoveGuidesWithFewReads(10);
            NGSResults.Log2Normalize_PerSample();

            NGSResults.AddData_Subtract(SampleBefore, SampleAfter, nDiff);
            NGSResults.AddData_Ranks(nDiff, nRankWI, nRankOA);

            Calculate_Individual_STARS();

            NGSResults.ExportGuides(NGSResults.OrigPath + "STARS.txt");
        }

        private void Calculate_Individual_STARS()
        {
            SortedList<double, List<gRNAClass>> sortedWithin;
            double STARS;
            foreach (GeneClass Gene in NGSResults.Genes.Values)
            {
                sortedWithin = new SortedList<double, List<gRNAClass>>(Gene.Guides.Count);
                foreach (gRNAClass Guide in Gene.Guides.Values)
                {
                    STARS = STARS_WorkOnGuide(Guide, Gene);
                    if (!sortedWithin.ContainsKey(STARS)) sortedWithin.Add(STARS, new List<gRNAClass>());
                    sortedWithin[STARS].Add(Guide);
                }
                for (int i = 0; i < sortedWithin.Values.Count; i++)
                    foreach (gRNAClass Guide in sortedWithin.Values[i])
                        Guide.AddData(nSTARSRank, i + 1);
            }
        }

        private double STARS_WorkOnGuide(gRNAClass guide, GeneClass Parent)
        {
            double STARS = Get_STARS_Value(guide, Parent, NGSResults.GuideCount);
            guide.AddData(nSTARS, STARS);
            return STARS;
        }

        public static double Get_STARS_Value(gRNAClass guide, GeneClass Parent, int TotalGuides)
        {
            double n = Parent.Guides.Count;
            double k = guide.GetData(nRankWI);
            double p = guide.GetData(nRankOA) / TotalGuides;
            double binomal = Math_Helper.BinomCoefficient((long)n, (long)k);
            double power1 = Math.Pow(p, k);
            double power2 = Math.Pow(1 - p, n - k);
            double result = binomal * power1 * power2;
            return result;
        }
    }

    public class gRNAClass
    {
        private Dictionary<string, double> _Data;
        public string Name { get; private set; }
        public string PostMod { get; private set; }
        public List<double> Reads { get; private set; }
        public double TotalReads { get; private set; }

        public gRNAClass(string Name, string PostMod)
        {
            Reads = new List<double>();
            _Data = new Dictionary<string, double>();
            TotalReads = 0;
            this.Name = Name;
            this.PostMod = PostMod;
        }

        public override string ToString()
        {
            return Name;
        }

        public void AddReads(double Reads)
        {
            this.Reads.Add(Reads);
            TotalReads += Reads;
        }

        public void AddData(string Key, double Value)
        {
            Key = Key.ToUpper();
            if (_Data.ContainsKey(Key))
                _Data[Key] = Value;
            else
                _Data.Add(Key, Value);
        }

        public double GetData(string Key)
        {
            Key = Key.ToUpper();
            if (_Data.ContainsKey(Key))
                return _Data[Key];
            else
                return double.NaN;
        }

        public string ExportLine()
        {
            return ExportLine(false, null);
        }

        public string ExportLine(bool Header, List<string> SampleNames)
        {
            StringBuilder sb = new StringBuilder();
            sb.Append(Header ? "Name" : Name); sb.Append("\t");
            sb.Append(Header ? "Mod" : PostMod); sb.Append("\t");
            for (int i = 0; i < Reads.Count; i++)
            {
                sb.Append(Header ? SampleNames[i] : Reads[i].ToString()); sb.Append("\t");
            }
            foreach (KeyValuePair<string, double> KVP in _Data)
            {
                sb.Append(Header ? KVP.Key : KVP.Value.ToString()); sb.Append("\t");
            }
            return sb.ToString();
        }
    }

    public class GeneClass
    {
        public string Name { get; private set; }
        public Dictionary<string, gRNAClass> Guides { get; private set; }

        public override string ToString()
        {
            return Name;
        }

        public GeneClass(string Name)
        {
            this.Name = Name;
            Guides = new Dictionary<string, gRNAClass>();
        }

        public gRNAClass AddOrGetGuide(string Name, string Mod)
        {
            string GK = GuideKey(Name, Mod);
            if (Guides.ContainsKey(GK))
                return Guides[GK];
            gRNAClass g = new gRNAClass(Name, Mod);
            AddGuide(g);
            return g;
        }

        public bool AddGuide(gRNAClass Guide)
        {
            string GK = GuideKey(Guide);
            if (Guides.ContainsKey(GK))
            {
                return false; //Already Has this One
            }
            Guides.Add(GK, Guide);
            return true;
        }

        public void RemoveGuideWithFewerThan(double MinimumReads)
        {
            List<string> ToRemove = new List<string>();
            foreach (gRNAClass g in Guides.Values)
            {
                if (g.TotalReads < MinimumReads)
                {
                    ToRemove.Add(GuideKey(g));
                }
            }
            for (int i = 0; i < ToRemove.Count; i++)
            {
                Guides.Remove(ToRemove[i]);
            }
        }

        public static string GuideKey(string Name, string Mod)
        {
            return Name + "_" + Mod;
        }

        public static string GuideKey(gRNAClass Guide)
        {
            return Guide.Name + "_" + Guide.PostMod;
        }
    }

    public class SampleClass
    {
        public string Name { get; internal set; }
        public double TotalReads { get; internal set; }

        public override string ToString()
        {
            return Name;
        }
        public SampleClass(string Name)
        {
            this.Name = Name;
        }

        public void AddReads(double Reads)
        {
            TotalReads += Reads;
        }
    }

    public static class Math_Helper
    {
        /// <summary>
        /// Calculates the binomial coefficient (nCk) (N items, choose k)
        /// </summary>
        /// <param name="n">the number items</param>
        /// <param name="k">the number to choose</param>
        /// <returns>the binomial coefficient</returns>
        public static long BinomCoefficient(long n, long k)
        {
            if (k > n) { return 0; }
            if (n == k) { return 1; } // only one way to chose when n == k
            if (k > n - k) { k = n - k; } // Everything is symmetric around n-k, so it is quicker to iterate over a smaller k than a larger one.
            long c = 1;
            for (long i = 1; i <= k; i++)
            {
                c *= n--;
                c /= i;
            }
            return c;
        }
    }
}
