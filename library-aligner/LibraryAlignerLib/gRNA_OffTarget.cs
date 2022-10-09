using System;
using System.Collections.Generic;
using System.Threading;
using System.Threading.Tasks;
using System.Text;
using System.IO;

namespace Library_Aligner
{
    namespace gRNA_OffTarget
    {

        public class Find_OffTargets
        {
            public static void A_StartHere()
            {
                Console.WriteLine("Generating MisMatch Index..");
                gRNA tGuide = new gRNA("CAGGACAGTAGTGATTCCTATGG");
                Find_OffTargets tFind = new Find_OffTargets(tGuide);

                Console.WriteLine("Loading Sequence(s)..");
                string Path = @"E:\f\Box Sync\Projects\GEiC\Customers\20171206 Mennerick\mm_ref_GRCm38.p4_chr11.fa";
                //string Path = @"E:\f\Box Sync\Projects\GEiC\Customers\20171206 Mennerick\rna.fa";
                tFind.LoadSequences(Path, false);

                Console.WriteLine("Getting other strand..");
                tFind.GenerateOtherStrand();

                Console.WriteLine("Finding OffTargets..");
                tFind.FindMatches();

                Console.WriteLine("Exporting..");
                tFind.ExportResults(Path + "_res.txt");
            }

            private Object thisLock = new object();
            private List<string> Results;
            private List<string> Lines;
            private List<string> LinesRev;
            private List<string> Headers;
            public int FlankingSequencesLength = 200;
            public gRNA Guide;

            public Find_OffTargets(gRNA GuideToSearch)
            {
                Guide = GuideToSearch;
            }

            public void ExportResults(string Path)
            {
                File.WriteAllLines(Path, Results);
            }

            /// <summary>
            /// Loads the sequences into the Find Class.
            /// </summary>
            /// <param name="FilePath">The Path of the FASTA File to Load</param>
            /// <param name="SeriesOfSequences">TRUE if this has many sequences separated by hearders >, false if this is just one long sequence (like a chromsome)</param>
            public void LoadSequences(string FilePath, bool SeriesOfSequences)
            {
                StreamReader SR = new StreamReader(FilePath);
                Lines = new List<string>();
                Headers = new List<string>();
                if (SeriesOfSequences)
                {
                    string more, h;
                    string t = SR.ReadLine();
                    while (!SR.EndOfStream)
                    {
                        if (t.StartsWith(">"))
                        {
                            h = t;
                            t = SR.ReadLine(); more = "";
                            while (!t.StartsWith(">"))
                            {
                                more += t;
                                if (SR.EndOfStream) break;
                                t = SR.ReadLine();
                            }
                            Lines.Add(more);
                            Headers.Add(h);
                        }
                    }
                }
                else
                {
                    Headers.Add(SR.ReadLine());
                    string t = SR.ReadToEnd();
                    Lines.Add(t.Replace("\n","").Replace("\r",""));
                }
            }

            public void GenerateOtherStrand()
            {
                LinesRev = new List<string>(Lines.Count);
                for (int i = 0; i < Lines.Count; i++)
                {
                    LinesRev.Add("");
                }
                Parallel.For(0, Lines.Count, (i) =>
                {
                    LinesRev[i] = SeqHelper.ReverseComplement(Lines[i]);
                });
            }

            public void FindMatches()
            {
                Results = new List<string>();
                AddGuide_Headers();
                bool UseParallel = Lines.Count > 48;
                //UseParallel = false; //ForDebugging
                if (UseParallel)
                    Parallel.For(0, Lines.Count, (i) =>
                     {
                         FindMatches_Single(Lines[i], LinesRev[i], Headers[i], !UseParallel);
                     });
                else
                    for (int i = 0; i < Lines.Count; i++)
                    {
                        //if (Headers[i].ToUpper().Contains("GABRG2")) {
                        FindMatches_Single(Lines[i], LinesRev[i], Headers[i], !UseParallel);
                    }
            }

            public void FindMatches_Single(string ForSeq, string RevSeq, string HeaderName, bool UseParallel)
            {
                int Max = ForSeq.Length - Guide.MinLength;
                //UseParallel = false; //ForDebugging
                if (UseParallel)
                    Parallel.For(0, Max, (position) =>
                    {
                        FindMatches_Single_Internal(ForSeq, RevSeq, position, Guide.MinLength, HeaderName);
                    });
                else
                    for (int position = 0; position < Max; position++)
                    {
                        FindMatches_Single_Internal(ForSeq, RevSeq, position, Guide.MinLength, HeaderName);
                    }
            }

            public bool FindMatches_Single_Internal(string ForSeq, string RevSeq, int Position, int MatchSize, string HeaderName)
            {
                bool found = false;
                string SubSeq = ForSeq.Substring(Position, MatchSize);
                for (int i = 0; i < 2; i++)
                {
                    if (i == 1) SubSeq = RevSeq.Substring(Position, MatchSize);
                    found = false;
                    if (Guide.Contains(SubSeq))
                    {
                        if (i == 1)
                        {

                        }
                        //If the Guide matches this smaller set, then look and see if there is a larger match as well (only add the larger)
                        if (MatchSize < Guide.MaxLength)
                        {
                            found = FindMatches_Single_Internal(ForSeq, RevSeq, Position, MatchSize + 1, HeaderName);
                            if (!found)
                                found = FindMatches_Single_Internal(ForSeq, RevSeq, i == 0 ? Position - 1 : Position + 1, MatchSize + 1, HeaderName);
                        }
                        if (!found)
                        {
                            Process_AddMatch(SubSeq, i == 0 ? ForSeq : RevSeq, Position, MatchSize, i == 0, HeaderName);
                        }
                        found = true;
                    }
                }
                return found;
            }

            private void Process_AddMatch(string SubSeq, string FullSeq, int Position, int MatchSize, bool SenseStrand, string HeaderName)
            {
                int Before_Start = Math.Max(Position - FlankingSequencesLength, 0);
                int Before_Len = Position - Before_Start;
                int After_Start = Position + MatchSize;
                int After_Len = Math.Min(After_Start + FlankingSequencesLength, FullSeq.Length) - After_Start;
                string Before = FullSeq.Substring(Before_Start, Before_Len);
                string After = FullSeq.Substring(After_Start, After_Len);
                AddGuide(SubSeq, Before, After, SenseStrand ? Position : (FullSeq.Length - Position), SenseStrand ? "+" : "-", HeaderName);
            }

            internal void AddGuide(string TargetSeq, string Before, string After, int Position, string Strand, string Name)
            {
                lock (thisLock)
                {
                    Results.Add(Name + "\t" + Position + "\t" + Strand + "\t" + Before + "\t" + TargetSeq + "\t" + After + "\t" + Guide.MainHash[TargetSeq].ToString());
                }
            }

            internal void AddGuide_Headers()
            {
                Results.Add("Name" + "\t" + "Pos" + "\t" + "Strand" + "\t" + "Before" + "\t" + "Target Seq" + "\t" + "After" + "\t" + "Score");
            }
        }

        public class gRNA
        {
            private static List<string> BPL = new List<string>(4) { "A", "C", "T", "G" };
            private Dictionary<int, List<List<byte>>> MisMatchPattern;
            private Dictionary<string, List<string>> MutationPattern;

            public Dictionary<string, byte> MainHash;
            public int MaxLength { get; internal set; }
            public int MaxBulge { get; internal set; }
            public int MaxMisMatches { get; internal set; }
            public int MinLength { get; internal set; }
            public string Sequence { get; internal set; }
            public bool Contains(string Seq) { return MainHash.ContainsKey(Seq); }


            public gRNA(string Sequence)
            {
                MaxBulge = 4; MaxMisMatches = 3;
                this.Sequence = Sequence;
                MaxLength = Sequence.Length;
                MinLength = MaxLength - MaxBulge + 1;
                Init_MMPattern();
                MutationPattern = new Dictionary<string, List<string>>();
                MainHash = new Dictionary<string, byte>(1500000);

                Generate_AllVariations(Sequence, MaxBulge, MaxMisMatches);
            }

            private void AddMutationPatterns(string SubSeq)
            {
                List<List<byte>> MMPattern = MisMatchPattern[SubSeq.Length];
                List<string> Seq = new List<string>(SubSeq.Length);
                for (int i = 0; i < SubSeq.Length; i++) Seq.Add(SubSeq.Substring(i, 1));
                HashSet<string> HS;
                List<string> tBP;
                foreach (List<byte> MMSet in MMPattern)
                {
                    HS = new HashSet<string>(); tBP = new List<string>();
                    for (int BP1 = 0; BP1 < 4; BP1++)
                    {
                        if (BPL[BP1] != Seq[MMSet[0]])
                        {
                            if (MMSet.Count > 1)
                            {
                                for (int BP2 = 0; BP2 < 4; BP2++)
                                {
                                    if (BPL[BP2] != Seq[MMSet[1]])
                                    {
                                        if (MMSet.Count > 2)
                                        {
                                            for (int BP3 = 0; BP3 < 4; BP3++)
                                            {
                                                if (BPL[BP3] != Seq[MMSet[2]])
                                                    MutPat_CheckAdd(SubSeq, MMSet, HS, BP1, BP2, BP3);
                                            }
                                        }
                                        else MutPat_CheckAdd(SubSeq, MMSet, HS, BP1, BP2, -1);
                                    }
                                }
                            }
                            else MutPat_CheckAdd(SubSeq, MMSet, HS, BP1, -1, -1);
                        }
                    }
                }
            }

            private void MutPat_CheckAdd(string Seq, List<byte> StartList, HashSet<string> Check, int BP1, int BP2, int BP3)
            {
                StringBuilder sb = new StringBuilder();
                if (BP1 > -1) sb.Append(BPL[BP1]);
                if (BP2 > -1) sb.Append(BPL[BP2]);
                if (BP3 > -1) sb.Append(BPL[BP3]);
                string t = sb.ToString();
                string key = MutPattern_Key(Seq, StartList);
                if (!Check.Contains(t))
                {
                    Check.Add(t);
                    if (!MutationPattern.ContainsKey(key))
                    {
                        MutationPattern.Add(key, new List<string>());
                    }
                    MutationPattern[key].Add(t);
                }
            }

            private string MutPattern_Key(string Seq, List<byte> MMSet)
            {
                if (MMSet.Count == 1) return Seq + ":" + MMSet[0].ToString();
                if (MMSet.Count == 2) return Seq + ":" + MMSet[0].ToString() + "_" + MMSet[1].ToString();
                if (MMSet.Count == 3) return Seq + ":" + MMSet[0].ToString() + "_" + MMSet[1].ToString() + "_" + MMSet[2].ToString();
                for (int i = 0; i < MMSet.Count; i++)
                {
                    //Not implemented
                }
                return "";
            }

            private string MutPat_Original(List<string> Seq, List<byte> MMSet)
            {
                StringBuilder sb = new StringBuilder();
                for (int i = 0; i < MMSet.Count; i++)
                {
                    sb.Append(Seq[MMSet[i]]);
                }
                return sb.ToString();
            }

            private void Init_MMPattern()
            {
                MisMatchPattern = new Dictionary<int, List<List<byte>>>();
                for (int i = MinLength; i <= MaxLength; i++)
                {
                    MisMatchPattern.Add(i, Init_MMPattern_Internal(i));
                }
            }

            private List<List<byte>> Init_MMPattern_Internal(int Len)
            {
                List<List<byte>> MMP = new List<List<byte>>();

                if (MaxMisMatches != 3)
                {
                    //NOT IMPLEMENTED
                }
                for (int i = 0; i < Len; i++)
                {
                    MMP.Add(new List<Byte>(1) { (byte)i });
                    for (int j = i + 1; j < Len; j++)
                    {
                        MMP.Add(new List<Byte>(2) { (byte)i, (byte)j });
                        for (int k = j + 1; k < Len; k++)
                        {
                            MMP.Add(new List<Byte>(3) { (byte)i, (byte)j, (byte)k });
                            /*for (int l = k + 1; l < Len - (MaxMisMatches - 4); l++) {
                                MisMatchPattern.Add(new List<Byte>(4) { (byte)i, (byte)j, (byte)k, (byte)l });
                            }*/
                        }
                    }
                }
                return MMP;
            }

            public void Generate_AllVariations(string Seq, int MaxBulge, int MaxMismatches)
            {
                AddMutationPatterns(Seq);
                Generate_Mismatches(Seq, MaxMismatches, 0); //First do it for the full sequence
                string Seq2;
                for (int Bulge_size = 1; Bulge_size < 4; Bulge_size++)
                {
                    //Parallel.For(0, (Seq.Length - Bulge_size - 2), Bulge_pos =>   //Tried this, but it actually made things worse
                    for (int Bulge_pos = 0; Bulge_pos < (Seq.Length - Bulge_size - 2); Bulge_pos++)
                    {
                        Seq2 = Seq.Substring(0, Bulge_pos) + Seq.Substring(Bulge_pos + Bulge_size);
                        AddMutationPatterns(Seq2);
                        Generate_Mismatches(Seq2, MaxMismatches, Bulge_size);
                    }
                }
            }

            public void Generate_Mismatches(string Seq, int Maxmismatches, int Bulge_size)
            {
                AddToHash(Seq, (byte)(Bulge_size * 10));
                char[] SeqChar = Seq.ToCharArray();
                string key;
                StringBuilder sb = new StringBuilder(Seq);
                List<List<byte>> MMPattern = MisMatchPattern[Seq.Length];
                List<string> MutsList;
                foreach (List<byte> MMSet in MMPattern)
                {
                    key = MutPattern_Key(Seq, MMSet);
                    MutsList = MutationPattern[key];
                    foreach (string Muts in MutsList)
                        AddToHash(GetSeq_MutationPattern(sb, SeqChar, MMSet, Muts), (byte)(Bulge_size * 10 + MMSet.Count));
                }
            }

            private string GetSeq_MutationPattern(StringBuilder sb, char[] Seq, List<byte> MMSet, string Mutations)
            {
                char[] MutChar = Mutations.ToCharArray();
                for (int i = 0; i < MMSet.Count; i++)
                {
                    sb[MMSet[i]] = MutChar[i];
                }
                string t = sb.ToString();
                //Reset it again
                for (int i = 0; i < MMSet.Count; i++)
                {
                    sb[MMSet[i]] = Seq[MMSet[i]];
                }
                return t;
            }

            public void Generate_Mismatches_Old(string Seq, int Maxmismatches, int Bulge_size)
            {
                AddToHash(Seq, (byte)(Bulge_size * 10));
                for (int Mismatches = 0; Mismatches < Maxmismatches; Mismatches++)
                {
                    AddToHash_AllVariations(Seq, Maxmismatches - Mismatches, Bulge_size);
                }
            }

            public void AddToHash_AllVariations(string Seq, int Mismatch_Depth, int Bulge_size)
            {
                string Seq2; string tBP;
                for (int Mismatch_pos = 0; Mismatch_pos < Seq.Length; Mismatch_pos++)
                {
                    for (int BP = 0; BP < 4; BP++)
                    {
                        tBP = BPL[BP];
                        if (Seq.Substring(Mismatch_pos, 1) != tBP)
                        {
                            Seq2 = Seq.Substring(0, Mismatch_pos) + tBP + Seq.Substring(Mismatch_pos + 1);
                            if (Mismatch_Depth > 1)
                                AddToHash_AllVariations(Seq2, Mismatch_Depth - 1, Bulge_size);
                            else
                                AddToHash(Seq2, (byte)((Bulge_size * 10) + 1));
                        }
                    }
                }
            }

            //private Object thisLock = new object();

            public bool AddToHash(string Seq, byte Score)
            {
                //lock (thisLock)  
                if (!MainHash.ContainsKey(Seq))
                {
                    //if (Seq.Length < CalcMinLength) CalcMinLength = Seq.Length;
                    MainHash.Add(Seq, Score);
                    return true;
                }
                return false;
            }
        }
    }
}
