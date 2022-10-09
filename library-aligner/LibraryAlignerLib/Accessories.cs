using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Library_Analysis;

namespace Library_Aligner
{
    public static class Accessories
    {
        private static string[] splitTab = new string[1] { "\t" };

        public static void Only_Join()
        {
            Console.WriteLine("Please wait... Joining");

            string PathRef = @"e:\temp\ngs\ngs\";

            Aligner A = new Aligner(null, null, "");
            Aligner.s.JoinRevType = 2;
            A.OnlyJoin(PathRef);
        }

        public static void Start_STARS_Analysis()
        {
            STARS_Analysis SA = new STARS_Analysis(@"E:\F\Box Sync\Projects\GEiC\CRISPR Library\20171025 Nishanth\Results\ref_w_res.txt", "0", "09-A05-", "10-A05-"); //"11-A05-"
            SA.RunAnalysis();
        }

        public static void Start_Protein_Aligner()
        {
            Console.WriteLine("Please wait...");
            string PathBase = @"E:\F";
            string PathRef = PathBase + @"\Box Sync\People\Jeff Milbrandt\Bacteria TIRs\20171029 PFAM\Raw\";
            string PathMatch = PathBase + @"\Box Sync\People\Jeff Gordon\92 Organisms\FAAs\";
            Library_Preference LP = new Library_Preference();
            Seq_List SL = new Seq_List();
            SL.Schema.AddUnique();
            SL.ImportAllFilesInFolder(PathRef, LP);
            Aligner A = new Aligner(SL, 5, PathRef, new Store());
            A.Preferences = LP;
            A.ExportSeqTable = true;
            A.AlignFolder(PathMatch, 7, 5, false);

            A.ExportResults(PathMatch + "res.txt");
            SL.Export_RefWithResults(A, PathMatch + "ref_w_res.txt", LP);
        }

        public static string Demux_gRNA_Sequences2Files()
        {
            System.Diagnostics.Debugger.Break();
            string PathRes = @"E:\F\Box Sync\Projects\GEiC\CRISPR Library\20180119 Asxl1\20180727\0727\res_f Poor Join.txt";
            PathRes = @"C:\Temp\NGS\NGS\res_1f.txt";
            FileInfo FI = new FileInfo(PathRes);
            DirectoryInfo DI = FI.Directory;
            DirectoryInfo DISub = DI.CreateSubdirectory("Subs");
            double ReadsToTargetMin = 28; //36;
            double ReadsToTarget;
            string Target; string t;
            int counter = 0;
            Dictionary<string, StreamWriter> WriterDict = new Dictionary<string, StreamWriter>();
            using (StreamReader SR = new StreamReader(PathRes))
            {
                string HeaderLine = SR.ReadLine();
                string[] arr;
                while (!SR.EndOfStream)
                {
                    t = SR.ReadLine();
                    arr = t.Split(splitTab, StringSplitOptions.None);
                    ReadsToTarget = double.Parse(arr[5]);
                    if (ReadsToTarget >= ReadsToTargetMin)
                    {
                        Target = arr[6];
                        if (!WriterDict.ContainsKey(Target))
                        {
                            WriterDict.Add(Target, new StreamWriter(Path.Combine(DISub.FullName, Target)));
                            WriterDict[Target].WriteLine(HeaderLine);
                        }
                        WriterDict[Target].WriteLine(t);
                        if (counter++ > 100000) { Console.WriteLine("Demuxing.."); counter = 0; }
                    }
                }
                SR.Close();
            }

            foreach (KeyValuePair<string, StreamWriter> KVP in WriterDict)
            {
                KVP.Value.Close();
            }

            return DISub.FullName;
        }

        public static void CountUMIsFolder(string FolderPath)
        {
            //Run Demux_gRNA_Sequences2Files and point this at the output folder
            //PathRes = @"E:\F\Box Sync\Projects\GEiC\CRISPR Library\20180119 Asxl1\20180727\0727\Subs\";
            StreamWriter SW = new StreamWriter(Path.Combine(FolderPath, "UMI_Out.txt"));
            DirectoryInfo DI = new DirectoryInfo(FolderPath);
            foreach (FileInfo FI in DI.GetFiles())
            {
                CountUMIs_File(FI.FullName, SW);
                Console.WriteLine(FI.Name);
            }
            SW.Close();
        }

        public static void CountUMIs_File(string FileName, StreamWriter OutPutFile)
        {
            Tuple<int, int> StartLen_5p, StartLen_3p;
            FindUMI_Locations(FileName, 1000, out StartLen_5p, out StartLen_3p);

            int Span = StartLen_3p.Item1 - StartLen_5p.Item1;
            if (Span < 85 && Span > 90) System.Diagnostics.Debugger.Break(); //Probably something went wrong
            if (StartLen_5p.Item2 != 5 || StartLen_3p.Item2 != 5) System.Diagnostics.Debugger.Break(); //Lengths aren't right

            Dictionary<string, Dictionary<string, double>> Track = new Dictionary<string, Dictionary<string, double>>();
            string t; string Headers; string[] arr = null; string UMI; double Occurences;
            using (StreamReader SR = new StreamReader(FileName))
            {
                Headers = SR.ReadLine(); //Headers
                while (!SR.EndOfStream)
                {
                    t = SR.ReadLine();
                    arr = t.Split(splitTab, StringSplitOptions.None);
                    UMI = GetUMI(arr[3], StartLen_5p, StartLen_3p);
                    Occurences = double.Parse(arr[1]);
                    if (!Track.ContainsKey(arr[0])) Track.Add(arr[0], new Dictionary<string, double>());
                    if (!Track[arr[0]].ContainsKey(UMI)) Track[arr[0]].Add(UMI, 0);
                    Track[arr[0]][UMI] += Occurences;
                }
                SR.Close();
            }

            foreach (KeyValuePair<string, Dictionary<string, double>> KVP in Track)
            {
                foreach (KeyValuePair<string, double> kUMI in KVP.Value)
                {
                    OutPutFile.WriteLine(FileName + "\t" + KVP.Key + "\t" + kUMI.Key + "\t" + kUMI.Value);
                }
            }
        }

        private static void FindUMI_Locations(string FileName, int MaxRead, out Tuple<int, int> StartLen_5p, out Tuple<int, int> StartLen_3p)
        {
            List<Dictionary<char, int>> Data = new List<Dictionary<char, int>>();
            string t;
            string[] arr = null;
            char[] carr; char c;
            Dictionary<char, int> PosDict;
            int ReadCount = 0;
            using (StreamReader SR = new StreamReader(FileName))
            {
                t = SR.ReadLine(); //Headers
                while (!SR.EndOfStream)
                {
                    t = SR.ReadLine();
                    arr = t.Split(splitTab, StringSplitOptions.None);
                    carr = arr[3].ToCharArray();
                    for (int i = 0; i < carr.Length; i++)
                    {
                        if (Data.Count <= i) Data.Add(new Dictionary<char, int>());
                        PosDict = Data[i];
                        c = carr[i];
                        if (!PosDict.ContainsKey(c)) PosDict.Add(c, 0);
                        PosDict[c]++;
                    }
                    if (ReadCount++ > MaxRead) break;
                }
                SR.Close();
            }

            //StringBuilder sB = new StringBuilder("i\tSumSq\r\n");
            double res; List<double> GinnyList = new List<double>(Data.Count);
            for (int i = 0; i < Data.Count; i++)
            {
                res = DictSumSq(Data[i]);
                GinnyList.Add(res);
                //sB.Append(i + "\t" + res + "\r\n");
            }
            //t = sB.ToString(); //Just for checking manually

            int GuidePos = int.Parse(arr[2]) + 10; //Middle of gRNA
            StartLen_5p = FindUMIsFrom_gRNA(GinnyList, GuidePos, -1);
            StartLen_3p = FindUMIsFrom_gRNA(GinnyList, GuidePos, +1);
        }

        public static string GetUMI(string AlignedSeq, Tuple<int, int> StartLen_5p, Tuple<int, int> StartLen_3p)
        {
            if (StartLen_3p.Item1 + StartLen_3p.Item2 > AlignedSeq.Length)
            {
                return "!ERROR";
            }
            return AlignedSeq.Substring(StartLen_5p.Item1, StartLen_5p.Item2) + "_" + AlignedSeq.Substring(StartLen_3p.Item1, StartLen_3p.Item2);
        }

        private static Tuple<int, int> FindUMIsFrom_gRNA(List<double> GinnyList, int Start, int Direction)
        {
            if (Direction < 0) Direction = -1; if (Direction >= 0) Direction = 1;
            int pointer = Start;
            double Stable = 0.9;
            double UMI = 0.3;
            int StretchStable = 1;
            int StretchUMI = 1;
            int Current = 0;
            int StartUMI = 0;
            int EndUMI = 0;
            while (true)
            {
                switch (Current) //Process
                {
                    case 0:
                        if (GinnyList[pointer] > Stable) StretchStable++;
                        if (GinnyList[pointer] < UMI) { Current = 1; StartUMI = pointer; }
                        break;
                    case 1:
                        if (GinnyList[pointer] > UMI) { EndUMI = pointer; pointer = Direction * GinnyList.Count; break; }
                        if (GinnyList[pointer] < UMI) StretchUMI++;
                        break;
                }
                //Move to Next
                pointer += Direction;
                if (Direction > 0) { if (pointer >= GinnyList.Count) break; }
                else { if (pointer < 0) break; }
            }
            if (Direction > 0)
                return new Tuple<int, int>(StartUMI, EndUMI - StartUMI);
            else
                return new Tuple<int, int>(EndUMI, StartUMI - EndUMI);
        }

        private static double DictSumSq(Dictionary<char, int> PosDict)
        {
            double d = 0; double sum = 0;
            foreach (int Val in PosDict.Values)
            {
                d += Val * Val;
                sum += Val;
            }
            return d / (sum * sum);
        }

        public static void Start_Demux()
        {
            Console.WriteLine("Please wait (Demux) ...");
            Library_Preference LP = new Library_Preference();
            Seq_List SL = new Seq_List(); SL.Schema.AddUnique();

            string PathRef = @"E:\Temp\NGS\MinIon\";
            string PathExport = @"E:\Temp\NGS\MinIon\Demux\";

            SL.ImportAllFilesInFolder(PathRef, LP);
            Aligner.s.MinSequencesForParallelAlignment = 10000;
            Aligner A;

            A = new Aligner(SL, 10, PathExport, null); //STORE needs to be added in this new version
            Aligner.s.AttemptJoin = false;
            Aligner.s.ExportLocation = false;
            Aligner.s.QualScoringSkip = true;
            A.Preferences = LP;
            A.DemuxFolder(PathRef, PathExport, 1);

        }

        public static void CombineFiles(string FolderPath)
        {
            FileInfo[] Files = new DirectoryInfo(FolderPath).GetFiles();
            bool first = true;
            StreamWriter sw = new StreamWriter(Path.Combine(FolderPath, "_Combined.txt"));
            string txt;
            foreach (FileInfo file in Files)
            {
                Console.WriteLine(file.Name);
                txt = File.ReadAllText(file.FullName);
                if (first)
                {
                    sw.Write(txt);
                    first = false;
                } else
                {
                    txt = txt.Substring(txt.IndexOf("\n") + 1);
                    sw.Write(txt);
                }
            }
            sw.Close();
        }
    }

    public class DistanceMatrix
    {
        private Dictionary<string, SeqWordRelationships> SeqRelHash;

        public DistanceMatrix()
        {
            SeqRelHash = new Dictionary<string, SeqWordRelationships>();
        }

        public static void aaaTry()
        {
            /*
            string s1 = "AGGAATGGTTTTTCAAAAACCTAATCCATTCCCTATGAGTATATATGATAATGTTGCATATGGACCAAGAATTCATGGTATAAAAGATAAAGCTACTTTAGACAAAATAGTTGAAGAGAGTTTAAGAGGAGCAGCTATTTGGGAAGAAGTAAAAGATAGATTAAATAAATCTGCACTAGGTTTATCTGGAGGGCAGCAACAACGTATATGTATAGCTAGAACTATAGCTATGAAGCCAGAAGTTATACTAATGGATGAGCCCACATCTGCATTAGATCCAATTTCTACATCTA";
            string s2 = "AGGAATGGTTTTTCAAAAACCTAATCCATTCCCTATGAGTATATATGATAATGTTGCATATGGACCAAGAATTCATGGTATAAAAGATAAAGCTACTTTAGATAAAATAGTTGAAGAAAGTTTAAGAGGAGCAGCTATTTGGGAAGAAGTAAAAGATAGATTAAATAAATCTGCATTAGGTTTATCTGGAGGGCAGCAACAACGTATATGTATAGCTAGAACTATAGCTATGAAGCCAGAAGTTATATTAATGGATGAGCCTACATCTGCACTAGACCCAATTTCTACATCTA";
            SeqWordRelationships swr1 = new SeqWordRelationships(s1, 2, 100);
            SeqWordRelationships swr2 = new SeqWordRelationships(s2, 2, 100);

            double inner;
            double outer;
            SeqDifference(swr1, swr2, out inner, out outer);
            */

            Console.WriteLine("Reading..");
            string Path = @"E:\F\Box Sync\People\Jeff Milbrandt\Bacteria TIRs\20171006 PPI\CompareS\Seqs.txt";
            //string Path = @"C:\zzrt\Box Sync\People\Jeff Milbrandt\Bacteria TIRs\20171006 PPI\CompareS\Seqs.txt";

            StreamReader SR = new StreamReader(Path);
            Dictionary<string, int> Seqs = new Dictionary<string, int>();
            string tSeq;
            while (!SR.EndOfStream)
            {
                tSeq = SR.ReadLine().ToUpper().Trim();
                if (Seqs.ContainsKey(tSeq))
                    Seqs[tSeq]++;
                else
                    Seqs.Add(tSeq, 1);
            }
            SR.Close();
            DistanceMatrix DM = new DistanceMatrix();
            List<string> DP2 = DM.Calculate_SeqDistances(Seqs.Keys.ToList<string>(), 2, 100);

            Console.WriteLine("Exporting..");
            Export(DP2, Path + "DP2.txt");

            /*
            List<List<double>> DM3;
            List<List<double>> DM4;
            DM3 = DM.CalculateDistanceMatrix(Seqs, 3);
            DM4 = DM.CalculateDistanceMatrix(Seqs, 4);

            Export(DM3, Seqs, @"E:\F\Box Sync\People\Jeff Milbrandt\Bacteria TIRs\20171006 PPI\CompareS\DM3.txt");
            Export(DM4, Seqs, @"E:\F\Box Sync\People\Jeff Milbrandt\Bacteria TIRs\20171006 PPI\CompareS\DM4.txt");
            */
        }

        public static void Export(List<List<double>> DistanceMatrix, List<string> Seq, string Path)
        {
            StreamWriter SW = new StreamWriter(Path);
            StringBuilder sb = new StringBuilder();
            sb.Append("\t");
            for (int i = 0; i < Seq.Count; i++)
            {
                sb.Append(Seq[i] + "\t");
            }
            SW.WriteLine(sb);
            for (int i = 0; i < DistanceMatrix.Count; i++)
            {
                sb = new StringBuilder();
                sb.Append(Seq[i]);
                for (int j = 0; j < DistanceMatrix[i].Count; j++)
                {
                    sb.Append(DistanceMatrix[i][j] + "\t");
                }
                SW.WriteLine(sb);
            }
            SW.Close();
        }

        public static void Export(List<string> PairList, string Path)
        {
            StreamWriter SW = new StreamWriter(Path);
            SW.WriteLine("Seq1\t" + "Seq2\t" + "Inner\t" + "Outer");
            for (int i = 0; i < PairList.Count; i++)
            {
                SW.WriteLine(PairList[i]);
            }
            SW.Close();
        }

        public List<string> Calculate_SeqDistances(List<string> Seqs, int Wordsize, int MaxDistance)
        {
            Console.WriteLine("Building SeqWordRelationships . . ");
            List<SeqWordRelationships> tList = new List<SeqWordRelationships>(Seqs.Count);
            Parallel.ForEach(Seqs, (Seq) =>
            {
                tList.Add(new SeqWordRelationships(Seq, Wordsize, MaxDistance));
            });

            Console.WriteLine("Calculating Differences . . ");
            List<string> Results = new List<string>(tList.Count);
            Seq_Heirachy sCluster = new Seq_Heirachy(tList.Count);

            Parallel.For(0, tList.Count - 1, i =>
            //for (int i = 0; i < tList.Count - 1; i++)
            {
                for (int j = i + 1; j < tList.Count; j++)
                {
                    double min; double max;
                    SeqDifference(tList[i], tList[j], out min, out max);
                    Results.Add(tList[i].Sequence + "\t" + tList[j].Sequence + "\t" + min + "\t" + max);
                    sCluster.AddResultPair(tList[i].Sequence, tList[j].Sequence, min, max);
                }
            });

            sCluster.Cluster_Sequences(0.04);
            return Results;
        }

        /*
        public List<List<double>> CalculateDistanceMatrix(List<string> Seqs, int Wordsize)
        {
            List<SeqWordRelationships> tList = new List<SeqWordRelationships>(Seqs.Count);
            Parallel.ForEach(Seqs, (Seq) => {
            //    foreach (string Seq in Seqs)
                    tList.Add(new SeqWordRelationships(Seq, Wordsize, 25));
            });

            //Prepare the list ahead of time
            List<List<double>> Results = new List<List<double>>(tList.Count);
            for (int i = 0; i < tList.Count - 1; i++)
            {
                Results.Add(new List<double>(tList.Count));
                for (int j = 0; j < tList.Count; j++) Results[i].Add(0);
            }

            double d = 0;
            Parallel.For(0, tList.Count - 1, i =>
            {
                for (int j = i + 1; j < tList.Count; j++)
                {
                    Results[i][j] = (SeqDifference(tList[i], tList[j]));
                }
            });

            return Results;
        }
        */

        public static double SeqDifference(string Seq1, string Seq2, int WordSize, int MaxDistance)
        {
            //Break up each sequence into pairs of words, measure the distances between the pairs
            SeqWordRelationships S1 = new SeqWordRelationships(Seq1, WordSize, MaxDistance);
            SeqWordRelationships S2 = new SeqWordRelationships(Seq2, WordSize, MaxDistance);
            double d; double i;
            SeqDifference(S1, S2, out i, out d);
            return d;
        }

        public static void SeqDifference(SeqWordRelationships SeqRels1, SeqWordRelationships SeqRels2, out double Inner, out double Outer)
        {
            if (SeqRels1.Sequence == SeqRels2.Sequence)
            {
                Inner = 1; Outer = 1; return;
            }
            SeqWordRelationships SLonger;
            SeqWordRelationships SShorter;
            if (SeqRels1.Length >= SeqRels2.Length)
            {
                SLonger = SeqRels1; SShorter = SeqRels2;
            }
            else
            {
                SLonger = SeqRels2; SShorter = SeqRels1;
            }
            //if (SLonger.Sequence.Contains(SShorter.Sequence)) { Inner = 1; Outer = (double)SShorter.Length / SLonger.Length; return; } //Not necessary I think

            //Get the list of common pairs from both sequences, also count total sets
            Inner = 0; Outer = 0;

            List<string> Pairs = new List<string>();
            foreach (string Pair in SLonger.InternalHash.Keys)
            {
                if (SShorter.InternalHash.ContainsKey(Pair))
                {
                    Pairs.Add(Pair);
                    //n += (SLonger.InternalHash[Pair].Count >= SShorter.InternalHash[Pair].Count) ? SLonger.InternalHash[Pair].Count : SShorter.InternalHash[Pair].Count;
                }
                else
                {
                    //n += SLonger.InternalHash[Pair].Count;
                }
            }
            //foreach (string Pair in SShorter.InternalHash.Keys) if (!SLonger.InternalHash.ContainsKey(Pair)) n += SShorter.InternalHash[Pair].Count;

            //Score each pairs similarity
            int good;
            int inner_short_bad;
            int outer_longer_bad;
            int sum_good = 0;
            int sum_is_bad = 0;
            int sum_ol_bad = 0;
            foreach (string Pair in Pairs)
            {
                CountCompare(SLonger.InternalHash[Pair], SShorter.InternalHash[Pair], out good, out inner_short_bad, out outer_longer_bad);
                sum_good += good;
                sum_is_bad += inner_short_bad;
                sum_ol_bad += outer_longer_bad;
            }
            Inner = (double)sum_good / SShorter.Elements;
            Outer = (double)sum_good / SLonger.Elements;
        }

        public static void CountCompare(Dictionary<int, int> PairDistancesLonger, Dictionary<int, int> PairDistancesShorter, out int Good, out int Inner_Bad, out int Outer_Bad)
        {
            Dictionary<int, int> PDL = new Dictionary<int, int>(PairDistancesLonger);
            Dictionary<int, int> PDS = new Dictionary<int, int>(PairDistancesShorter);
            Good = 0; Inner_Bad = 0;
            Outer_Bad = 0;
            foreach (KeyValuePair<int, int> KVP in PairDistancesShorter)
            {
                if (CountCompare_Internal(PairDistancesLonger, KVP.Key, KVP.Value, 1, PDS, PDL, ref Good, ref Inner_Bad, ref Outer_Bad))
                {
                    //All taken care of
                }
                else
                {
                    //Just need to turn this up a bit
                    Inner_Bad += KVP.Value;
                }
            }
            /*
            bool found;
            foreach (KeyValuePair<int, int> KVP in PDS)
            {
                found = false;
                if (CountCompare_Internal(PDL, KVP.Key + 1, KVP.Value, 0.5, PDS, PDL, ref Good, ref Inner_Bad, ref Outer_Bad)) found = true;
                if (!found)
                    if (CountCompare_Internal(PDL, KVP.Key - 1, KVP.Value, 0.5, PDS, PDL, ref Good, ref Inner_Bad, ref Outer_Bad)) found = true;
                if (found)
                {

                }
            } 
            */
            foreach (KeyValuePair<int, int> KVP in PDL)
            {
                Outer_Bad += KVP.Value;
            }
        }


        private static bool CountCompare_Internal(Dictionary<int, int> PairDistancesLonger, int Key, int Value, double weight, Dictionary<int, int> PDS, Dictionary<int, int> PDL, ref int Good, ref int Inner_Bad, ref int Outer_Bad)
        {
            int tGood;
            int i; int o;
            if (PairDistancesLonger.ContainsKey(Key))
            {
                i = Value;
                o = PairDistancesLonger[Key];
                tGood = (int)(Math.Min(i, o) * weight);
                Good += tGood;
                Inner_Bad += Math.Max(0, i - o);
                Outer_Bad += Math.Max(0, o - i);
                PDS[Key] -= tGood; if (PDS[Key] == 0) PDS.Remove(Key);
                PDL[Key] -= tGood; if (PDL[Key] == 0) PDL.Remove(Key);
                return true;
            }
            return false;
        }

        public static double FractionCompare(List<int> PairDistances1, List<int> PairDistances2)
        {
            List<int> PairsMore;
            List<int> PairsLess;
            if (PairDistances1.Count >= PairDistances2.Count)
            {
                PairsMore = PairDistances1; PairsLess = PairDistances2;
            }
            else
            {
                PairsMore = PairDistances2; PairsLess = PairDistances1;
            }
            double d = 0; double best = 0; double sum = 0;
            for (int i = 0; i < PairsMore.Count; i++)
            {
                best = double.MaxValue;
                for (int j = 0; j < PairsLess.Count; j++)
                {
                    d = Math.Abs(PairsMore[i] - PairsLess[j]) + 1;
                    if (d < best) best = d;
                }
                sum += 1 / (best);
            }
            return sum;
        }
    }

    public class Seq_Heirachy
    {
        public Dictionary<string, double> BigDistanceHash;
        public Dictionary<string, DistanceMatrix_ResultsSet> BigDMRPHash;
        public List<DistanceMatrix_ResultsSet> DMResPairs;
        public Dictionary<string, List<DistanceMatrix_ResultsSet>> PairHash;
        public Dictionary<string, ClusterNode> Nodes;
        public Object thisLock = new object();

        private Dictionary<string, int> _SequenceOrderHash;
        private int _SequenceOrderTrack;

        public Seq_Heirachy(int Capacity)
        {
            DMResPairs = new List<DistanceMatrix_ResultsSet>(Capacity);
            BigDistanceHash = new Dictionary<string, double>(Capacity);
            BigDMRPHash = new Dictionary<string, DistanceMatrix_ResultsSet>(Capacity);
            PairHash = new Dictionary<string, List<DistanceMatrix_ResultsSet>>();
            _SequenceOrderHash = new Dictionary<string, int>();
            _SequenceOrderTrack = 0;
        }

        public void AddResultPair(string Seq1, string Seq2, double inner_min, double outer_max)
        {
            DistanceMatrix_ResultsSet DMRP = new DistanceMatrix_ResultsSet(Seq1, Seq2, inner_min, outer_max);
            DMResPairs.Add(DMRP);
            lock (thisLock)
            {
                foreach (string Seq in DMRP.SeqList)
                {
                    if (!PairHash.ContainsKey(Seq))
                        PairHash.Add(Seq, new List<DistanceMatrix_ResultsSet>());
                    PairHash[Seq].Add(DMRP);

                    if (!_SequenceOrderHash.ContainsKey(Seq)) _SequenceOrderHash.Add(Seq, _SequenceOrderTrack++);
                }
                string BDHK = BigDistHashKey(DMRP.SeqList[0], DMRP.SeqList[1]);
                if (!BigDistanceHash.ContainsKey(BDHK)) BigDistanceHash.Add(BDHK, DMRP.Distance);
                if (!BigDMRPHash.ContainsKey(BDHK)) BigDMRPHash.Add(BDHK, DMRP);
            }
        }

        public string BigDistHashKey(string Seq1, string Seq2)
        {
            if (_SequenceOrderHash[Seq1] < _SequenceOrderHash[Seq2])
                return Seq1 + "_" + Seq2;
            else
                return Seq2 + "_" + Seq1;
        }

        public void Cluster_Sequences(double MaxDistanceSiblingVsChild)
        {
            //First get out the best sequence pair
            List<DistanceMatrix_ResultsSet> Best = new List<DistanceMatrix_ResultsSet>();
            Best = LowestDistance(DMResPairs);
            DistanceMatrix_ResultsSet tDMRP = Best[0];
            ClusterNode CN1 = ClusterNode.NewBranch(tDMRP);

            //Figure out if there is anything else

            for (int i = 0; i < tDMRP.SeqList.Count; i++)
            {
                foreach (DistanceMatrix_ResultsSet DMRP in PairHash[tDMRP.SeqList[i]])
                {
                    if ((DMRP != tDMRP) && (!CN1.ContainsAllSequences(DMRP)) && (DMRP.Distance < MaxDistanceSiblingVsChild))
                    {
                        //Ok, we know this one's distance is low, what about its relationship to the new member?
                        string SeqNew = tDMRP.SeqNotInCommon(DMRP);
                        string SeqOrig = DMRP.SeqNotInCommon(tDMRP);
                        string BDHK = BigDistHashKey(SeqNew, SeqOrig);
                        double distanceNew = BigDistanceHash[BDHK];
                        if (distanceNew < MaxDistanceSiblingVsChild)
                            CN1.AddChild(SeqNew, distanceNew, DMRP.Distance);
                    }
                }
            }

            //Re-build the BigDistance Matrix with this information
            foreach (KeyValuePair<string, int> Seq in CN1.SeqHash)
            {
                foreach (DistanceMatrix_ResultsSet DMRP in PairHash[Seq.Key])
                {
                    string SeqNew = CN1.SeqNotInCommon(DMRP);
                    if (SeqNew != string.Empty)
                    {
                        List<double> dists = new List<double>(CN1.SeqHash.Count);
                        double SumSq = 0;
                        foreach (KeyValuePair<string, int> SeqInternal in CN1.SeqHash)
                        {
                            string BDHK = BigDistHashKey(SeqNew, SeqInternal.Key);
                            if (BigDistanceHash.ContainsKey(BDHK))
                            {
                                dists.Add(BigDistanceHash[BDHK]);
                                SumSq += Math.Pow(BigDistanceHash[BDHK], 2);
                                BigDistanceHash.Remove(BDHK);
                            }
                        }
                        SumSq = Math.Sqrt(SumSq);
                        BigDistanceHash.Add("NewSeq", SumSq);
                    }
                }
            }
        }

        public List<DistanceMatrix_ResultsSet> LowestDistance(List<DistanceMatrix_ResultsSet> Check)
        {
            List<DistanceMatrix_ResultsSet> Best = new List<DistanceMatrix_ResultsSet>();
            double MinDist = double.MaxValue;
            foreach (DistanceMatrix_ResultsSet DMRP in Check)
            {
                if (DMRP.Distance < MinDist)
                {
                    MinDist = DMRP.Distance;
                    Best = new List<DistanceMatrix_ResultsSet>();
                }
                if (DMRP.Distance == MinDist)
                    Best.Add(DMRP);
            }
            return Best;
        }

        public void Cluster_Sequences_PoorMethod(double Param1)
        {
            //SortedList<double, List<DistanceMatrix_ResultsPair>> sorted = new SortedList<double, List<DistanceMatrix_ResultsPair>>();
            while (DMResPairs.Count > 4)
            {
                List<DistanceMatrix_ResultsSet> Best = new List<DistanceMatrix_ResultsSet>();
                Best = LowestDistance(DMResPairs);

                DistanceMatrix_ResultsSet tDMRP = Best[0];
                DMResPairs.Remove(tDMRP);
                int count = 0;
                foreach (DistanceMatrix_ResultsSet DMRP in DMResPairs)
                {
                    if (DMRP.Contains(tDMRP))
                    {
                        DMRP.ReStructure(tDMRP, Param1);
                        count++;
                    }
                }
            }
        }
    }

    public class ClusterNode
    {
        public List<double> Distances;
        public Dictionary<string, int> SeqHash;
        public List<ClusterNode> ChildNodes;
        //public Dictionary<ClusterNode, int> ChildNodes;
        public bool Leaf;

        private void Init()
        {
            Distances = new List<double>();
            SeqHash = new Dictionary<string, int>();
            ChildNodes = new List<ClusterNode>();
        }

        public ClusterNode()
        {
            Init();
        }

        public ClusterNode(DistanceMatrix_ResultsSet DMRP, int Index)
        {
            Init();
            Leaf = true;
            if (SeqHash.ContainsKey(DMRP.SeqList[Index]))
            {

            }
            else
            {
                SeqHash.Add(DMRP.SeqList[Index], 0);
            }
        }

        public ClusterNode AddChild(ClusterNode Child)
        {
            Leaf = false;
            ChildNodes.Add(Child);
            AddChildSequences(Child);
            return ChildNodes[ChildNodes.Count - 1];
        }

        public ClusterNode AddChild(string NewSequence, double NewDistance1, double NewDistance2)
        {
            ClusterNode CN = new ClusterNode();
            CN.SeqHash.Add(NewSequence, 0);
            AddChild(CN);
            Distances.Add(NewDistance1);
            Distances.Add(NewDistance2);
            return CN;
        }

        private void AddChildSequences(ClusterNode Child)
        {
            foreach (KeyValuePair<string, int> KVP in Child.SeqHash)
            {
                SeqHash.Add(KVP.Key, KVP.Value + 1);
            }
        }

        public bool ContainsAnySequences(DistanceMatrix_ResultsSet DMRP)
        {
            foreach (string Seq in DMRP.SeqList)
            {
                if (SeqHash.ContainsKey(Seq)) return true;
            }
            return false;
        }

        public bool ContainsAllSequences(DistanceMatrix_ResultsSet DMRP)
        {
            foreach (string Seq in DMRP.SeqList)
            {
                if (!SeqHash.ContainsKey(Seq)) return false;
            }
            return true;
        }

        public string SeqNotInCommon(DistanceMatrix_ResultsSet DMRP)
        {
            foreach (string Seq in DMRP.SeqList)
            {
                if (!SeqHash.ContainsKey(Seq)) return Seq;
            }
            return string.Empty;
        }

        public static ClusterNode NewBranch(DistanceMatrix_ResultsSet DMRP)
        {
            ClusterNode tCN = new ClusterNode();
            tCN.Leaf = false;
            tCN.AddChild(new ClusterNode(DMRP, 0));
            tCN.AddChild(new ClusterNode(DMRP, 1));
            tCN.Distances.Add(DMRP.Distance);
            return tCN;
        }
    }

    public class DistanceMatrix_ResultsSet
    {
        internal double _Inner;
        internal double _Outer;
        internal double _Distance;
        internal bool _Dirty;

        internal List<DistanceMatrix_ResultsSet> _ChildSet;
        internal List<string> _Seq;
        internal Dictionary<string, bool> _SeqDict;
        internal List<double> _Distances;

        public List<double> Distances { get => _Distances; }

        public double Inner
        {
            get
            {
                Refresh();
                return _Inner;
            }
        }

        public double Outer
        {
            get
            {
                Refresh();
                return _Outer;
            }
        }

        public double Distance
        {
            get
            {
                Refresh();
                return _Distance;
            }
        }

        public List<string> SeqList { get => _Seq; }

        public List<DistanceMatrix_ResultsSet> ChildSet { get => _ChildSet; }

        public DistanceMatrix_ResultsSet this[int index]
        {
            get
            {
                return ChildSet[index];
            }
        }

        private void Init()
        {
            _Seq = new List<string>();
            _SeqDict = new Dictionary<string, bool>();
            _ChildSet = new List<DistanceMatrix_ResultsSet>();
            _Distances = new List<double>();
        }

        public DistanceMatrix_ResultsSet(string seq1, string seq2, double inner, double outer)
        {
            Init();

            AddSequence(seq1, true);
            AddSequence(seq2, true);

            _Inner = inner;
            _Outer = outer;
            _Distance = 1 - outer;
            _Distances.Add(_Distance);
            _Dirty = false;
        }

        /// <summary>
        /// Adds a new sequence to this Results Set
        /// </summary>
        /// <param name="Seq">String of the Sequence</param>
        /// <param name="IsOriginal">This is TRUE if this is the leaf node that contains this sequence, false if adding as part of another results set</param>
        public void AddSequence(string Seq, bool IsOriginal)
        {
            //Only adds if it is not already here
            if (!_SeqDict.ContainsKey(Seq))
            {
                _SeqDict.Add(Seq, IsOriginal);
                _Seq.Add(Seq);
            }
        }

        public void Refresh()
        {
            if (_Dirty)
            {
                double d = 0;
                foreach (double dist in _Distances)
                {
                    d += dist;
                }
                _Dirty = false;
            }
        }

        public bool Contains(DistanceMatrix_ResultsSet DMRP)
        {
            for (int i = 0; i < DMRP.SeqList.Count; i++)
            {
                if (Contains(DMRP.SeqList[i])) return true;
            }
            return false;
        }

        public bool Contains(string Sequence)
        {
            return _SeqDict.ContainsKey(Sequence);
        }

        public string SeqInCommon(DistanceMatrix_ResultsSet DMRP)
        {
            for (int i = 0; i < DMRP.SeqList.Count; i++)
            {
                if (Contains(DMRP.SeqList[i])) return DMRP.SeqList[i];
            }
            return string.Empty;
        }

        public string SeqNotInCommon(DistanceMatrix_ResultsSet DMRP)
        {
            for (int i = 0; i < DMRP.SeqList.Count; i++)
            {
                if (!Contains(DMRP.SeqList[i])) return DMRP.SeqList[i];
            }
            return string.Empty;
        }

        public void ReStructure(DistanceMatrix_ResultsSet DMRP, double SiblingVsChild)
        {
            string SeqMatch = SeqInCommon(DMRP);
            if (SeqMatch != string.Empty)
            {
                if (Distance <= SiblingVsChild)
                {
                    //if the distance is low enough, then we will want to add this as a sibling rather than a child
                    if (_SeqDict[SeqMatch])
                    {

                    }
                    else
                    {

                    }
                }
                if (_SeqDict[SeqMatch])
                {
                    //Replace whatever sequences matches DMRP with this
                    _ChildSet.Add(DMRP);
                }
                else
                {
                    //Add this to a child rather than to me
                    foreach (DistanceMatrix_ResultsSet Child in _ChildSet)
                    {
                        if (Child.Contains(DMRP))
                        {
                            Child.ReStructure(DMRP, SiblingVsChild);
                            break;
                        }
                    }
                }
                _Distances.AddRange(DMRP.Distances);
                foreach (string Seq in DMRP.SeqList)
                {
                    if (Contains(Seq))
                        _SeqDict[Seq] = false; //This is no longer the unique one
                    else
                        AddSequence(Seq, false);
                }
                _Dirty = true;
            }
        }

    }

    public class SeqWordRelationships
    {
        public Dictionary<string, Dictionary<int, int>> InternalHash;
        public string Sequence;
        public int WordSize;
        public int MaxDistance;
        public int Length;
        public int Elements;

        public SeqWordRelationships(string Sequence)
        {
            Init(Sequence, 3, 25);
        }

        public SeqWordRelationships(string Sequence, int WordSize)
        {
            Init(Sequence, WordSize, 25);
        }

        public SeqWordRelationships(string Sequence, int WordSize, int MaxDistance)
        {
            Init(Sequence, WordSize, MaxDistance);
        }

        private void Init(string Sequence, int WordSize, int MaxDistance)
        {
            this.Sequence = Sequence;
            this.WordSize = WordSize;
            this.MaxDistance = MaxDistance;
            this.Length = Sequence.Length;

            BreakSeq(Sequence, WordSize, MaxDistance, out InternalHash, out Elements);
        }

        public static void BreakSeq(string Seq, int WordSize, int MaxDistance, out Dictionary<string, Dictionary<int, int>> tHash, out int Elements)
        {
            tHash = new Dictionary<string, Dictionary<int, int>>();
            List<string> tWords = new List<string>(Seq.Length - WordSize);
            string t = "";
            Elements = 0;
            for (int i = 0; i < Seq.Length - WordSize + 1; i++)
            {
                t = Seq.Substring(i, WordSize);
                tWords.Add(t);
            }
            int dist = 0;
            for (int i = 0; i < tWords.Count - 1; i++)
            {
                for (int j = i + 1; j < tWords.Count; j++)
                {
                    dist = j - i;
                    if (dist <= MaxDistance)
                    {
                        t = tWords[i] + "_" + tWords[j];
                        if (!tHash.ContainsKey(t))
                            tHash.Add(t, new Dictionary<int, int>());
                        if (!tHash[t].ContainsKey(dist))
                            tHash[t].Add(dist, 0);
                        tHash[t][dist]++;
                        Elements++;
                    }
                }
            }
        }
    }

    public static class SeqHelper
    {
        public static Dictionary<char, char> HashComplement = new Dictionary<char, char>(4) { { 'C', 'G' }, { 'G', 'C' }, { 'A', 'T' }, { 'T', 'A' } };

        public static string ReverseComplement(string Sequence)
        {
            StringBuilder sb = new StringBuilder();
            char c;
            char[] Seq = Sequence.ToUpper().Trim().ToCharArray();
            for (int i = Seq.Length - 1; i >= 0; i--)
            {
                c = Seq[i];
                if (HashComplement.ContainsKey(c))
                    sb.Append(HashComplement[c]);
                else
                    sb.Append(c);
            }
            return sb.ToString();
        }
    }

    public class AlignedResultsManager
    {
        public StreamWriter PrimaryResults;
        public ConcurrentDictionary<StreamWriter, object> DemuxLocking;
        public ConcurrentDictionary<string, StreamWriter> DemuxResults;
        public string FolderPath;
        public string FolderPathDemux;
        public static Dictionary<string, StreamWriter> OpenStreams = new Dictionary<string, StreamWriter>();

        public AlignedResultsManager(string FolderPath, Seq_List SL)
        {
            this.FolderPath = FolderPath;
            FileInfo FI = new FileInfo(Path.Combine(FolderPath, "res_1f.txt"));
            if (OpenStreams.ContainsKey(FI.FullName))
            {
                PrimaryResults = OpenStreams[FI.FullName];
                if (PrimaryResults.BaseStream ==  null)
                {
                    //This means it was already closed, shouldn't be here . .
                }
            }
            else
            {
                PrimaryResults = new StreamWriter(FI.FullName /*, FI.Exists*/);
                OpenStreams.Add(FI.FullName, PrimaryResults);
                ExportPrimary(null); //Headers
            }
            if (Aligner.s.Demux)
            {
                this.FolderPathDemux = Path.Combine(FolderPath, "Demux");
                DirectoryInfo DI = new DirectoryInfo(this.FolderPathDemux);
                if (!DI.Exists) DI.Create();
                DemuxResults = new ConcurrentDictionary<string, StreamWriter>();
                DemuxLocking = new ConcurrentDictionary<StreamWriter, object>();
                foreach (Seq seq in SL.List)
                {
                    StreamWriter SW = new StreamWriter(Path.Combine(FolderPathDemux, seq.Name.Replace("|","") + ".txt"));
                    DemuxResults.TryAdd(seq.Name, SW);
                    DemuxLocking.TryAdd(SW, new object());
                    ExportDemux(SW); //Headers
                }
            }
        }

        public void ExportResults(AlignResultsClass ARC)
        {
            ExportPrimary(ARC);
            if (Aligner.s.Demux)
            {
                if (ARC.Target != "")
                {
                    if (!DemuxResults.ContainsKey(ARC.Target))
                    {
                        System.Diagnostics.Debugger.Break(); //These should all be built ahead of time, so something went wrong :(
                                                             //if (ARC.Target.Length < 5) System.Diagnostics.Debugger.Break();
                                                             //StreamWriter SW = new StreamWriter(Path.Combine(FolderPathDemux, ARC.Target + ".txt")); //This only works if NOT Parallel
                                                             //DemuxResults.GetOrAdd(ARC.Target, new StreamWriter(Path.Combine(FolderPathDemux, ARC.Target + ".txt")));
                                                             //DemuxLocking.TryAdd(SW, new object());
                    }
                    ExportDemux(ARC);
                }
            }
        }

        private static object PrimaryLock = new object();

        public void ExportPrimary(AlignResultsClass ARC)
        {
            lock (PrimaryLock)
            {
                if (ARC == null)
                {
                    PrimaryResults.WriteLine("Pre" + "\t" + "Prefix" + "\t" + "Occurences" + "\t" + "Good Qual Length" + "\t" + "Common Matches" + "\t" + "Reads to Target" + "\t" + "Target" + "\t" + "Other Target");
                }
                else
                {
                    PrimaryResults.WriteLine("" + "\t" + ARC.Prefix + "\t" + ARC.occurences + "\t" + ARC.seqLength + "\t" + ARC.RepresentsAllCount + "\t" + ARC.Max + "\t" + ARC.Target + "\t" + "" + "\t");
                }
            }
        }

        public void ExportDemux(StreamWriter SW)
        {
            ExportDemux(SW, null);
        }

        public void ExportDemux(AlignResultsClass ARC)
        {
            ExportDemux(null, ARC);
        }

        public void ExportDemux(StreamWriter SW, AlignResultsClass ARC)
        {
            if (SW == null) SW = DemuxResults[ARC.Target];
            object tLock = DemuxLocking[SW];
            lock (tLock)
            {
                if (ARC == null)
                {
                    SW.WriteLine("Prefix" + "\t" + "Occurences" + "\t" + "Start" + "\t" + "AlignedSequence");
                }
                else
                {
                    if (ARC.Max >= 5)
                    {
                        SW.WriteLine(ARC.Prefix + "\t" + ARC.occurences + "\t" + ARC.PositionInAligned + "\t" + ARC.AlignedSequence());
                    }
                }
            }
        }

        public void CloseFiles()
        {
            //20210601 Update below, before just closed the primary stream
            foreach (StreamWriter SR in OpenStreams.Values) SR.Close();
            OpenStreams = new Dictionary<string, StreamWriter>();
            
            if (Aligner.s.Demux)
            {
                foreach (StreamWriter SW in DemuxResults.Values)
                {
                    SW.Close();
                }
            }
        }
    }

    public class AlignResultsClass
    {
        public string Prefix;
        public double occurences;
        public string seq;
        public int seqLength;
        public int RepresentsAllCount;
        public double Max;
        public int Location;
        public List<string> TargetList;

        public static int Global_SpaceBeforeTargetSeq = -1;
        public static int Global_Buffer = 60;

        public string Target
        {
            get
            {
                return TargetList.Count > 0 ? TargetList[0] : "";
            }
        }

        public int PositionInAligned
        {
            get
            {
                if (Global_SpaceBeforeTargetSeq == -1)
                {
                    //First time we are doing this so we have to chose the total distance
                    Global_SpaceBeforeTargetSeq = Global_Buffer + this.Location;
                }
                return Global_SpaceBeforeTargetSeq;
            }
        }

        public AlignResultsClass(string Prefix, double occurences, string seq, int seqLength, int RepresentsAllCount, double Max)
        {
            TargetList = new List<string>();
            this.Prefix = Prefix;
            this.occurences = occurences;
            this.seq = seq;
            this.seqLength = seqLength;
            this.RepresentsAllCount = RepresentsAllCount;
            this.Max = Max;

        }

        public string ExportOld()
        {
            return Prefix + "\t" + occurences + "\t" + seq + "\t" + seqLength + "\t" + RepresentsAllCount + "\t" + Max + "\t" + Target;
        }

        public string AlignedSequence()
        {
            if (Global_SpaceBeforeTargetSeq == -1)
            {
                //First time we are doing this so we have to chose the total distance
                Global_SpaceBeforeTargetSeq = Global_Buffer + this.Location;
            }
            int adjust = Global_SpaceBeforeTargetSeq - this.Location;
            if (adjust < 0)
            {
                return this.seq.Substring(-adjust);
            }
            return new string(' ', adjust) + this.seq;
        }
    }

}
