using System.IO;
using System.Linq;
using System.Collections.Generic;
using System.Xml.Serialization;
using System.Collections.Concurrent;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using System;
using System.Collections;

namespace Library_Aligner
{

    public static class Library_Aligner
    {
        private static bool _ConsoleExport = false;

        public static void Start_Library_Aligner(string XMLSettingsPath = "", bool ConsoleExport = false)
        {
            Store SToUse;
            //Check to see if there is a settings file, if not use detaults
            if (XMLSettingsPath == "")
            {
                Message("Using Defaults since no settings file was specified.", false, false);
                SToUse = Store.GetDefaults();
            }
            else SToUse = Load(XMLSettingsPath);
            Start_Library_Aligner(SToUse, null, ConsoleExport);
        }

        public static void Start_Library_Aligner(Store StoreToUse, System.ComponentModel.BackgroundWorker BW = null, bool ConsoleExport = false)
        {
            S = StoreToUse;
            Seq_List.SkipRankings = true;
            _ConsoleExport = ConsoleExport;

            if (ConsoleExport)
            {
                Console.WriteLine("");
                Console.WriteLine(" - - Buchser Lab / GEiC 'Library Aligner' 2017 - - ");
                Console.WriteLine("1: Reference sequences (tab delimited) name (tab) sequence. 2: FastQ sequences.  3: (optional) Settings file.");
                Console.WriteLine("");
                Console.WriteLine("Please wait...");
            }
            //Setup Some of the Objects
            Library_Preference LP = new Library_Preference(); Aligner A = null; Seq_List SL = new Seq_List(); SL.Schema.AddUnique(); SetupSchemas();

            //Check to see if the paths Exist
            if (!new DirectoryInfo(S.PathFastQToSearch).Exists) { Message("Couldn't find the folder where the sequences are supposed to be.", true, false, BW); return; }
            if (!new DirectoryInfo(S.PathReferenceList).Exists) { Message("Couldn't find Sequences to Search for.  Make a Tab-Delimited Text file with the first column the identified and the second column the sequence.", true, true, BW); return; }
            string Note = SL.ImportAllFilesInFolder(S.PathReferenceList, LP); if (Note != "") { Message(Note, true, false, BW); return; }

            //Import the sequences that you are searching for
            if (SL.List.Count < 1) { Message("Didn't find any sequences. Please check the ref path provided.", true, false, BW); return; }
            //if (SL.List.Count < 3) { Message("Only found " + SL.List.Count + " sequences . . if this is correct press a key to continue, if not close the program.", true, true); }
            if (BW != null) BW.ReportProgress(0, "Found " + SL.List.Count.ToString("0,000") + " reference sequences. Libraries = " + string.Join(", ", SL.Libraries));

            //Check the AlignSize
            FinalAlignSize = S.AlignSize;
            if (FinalAlignSize < 0) FinalAlignSize = SL.MeanAlignSize;
            if (FinalAlignSize < S.PieceSize) FinalAlignSize = S.PieceSize;
            //if (FinalAlignSize > 150) FinalAlignSize = 150;

            //Only for advanced
            bool RunOptimization = false; if (RunOptimization) A = RunOptimizations(A, SL);

            //Here is where the Aligner actually starts
            A = new Aligner(SL, S.PieceSize, S.PathFastQToSearch, S); A.Preferences = LP;
            if (S.SplitLargeFiles) A.SplitLargeFiles(S.PathFastQToSearch, 2000);

            //Align the folder
            A.AlignFolder(S.PathFastQToSearch, MinReads2Target, S.MinutesMaxPerIndex, true, BW);

            //Get ready for Export
            if (!Seq_List.SkipRankings) SL.CalculateRanks(A, LP);
            if (ConsoleExport) { Console.WriteLine(); Console.WriteLine("Exporting . . "); } else { if (BW != null) BW.ReportProgress(80, "Exporting . . "); }
            A.ExportResults(Path.Combine(S.PathFastQToSearch, "res_f.txt"));
            SL.Export_RefWithResults(A, Path.Combine(S.PathFastQToSearch, "ref_f.txt"), LP); //Sometimes Spotfire thinks that important columns are integer . . make sure to set them as reals
            SL.Export_RefLong(A, Path.Combine(S.PathFastQToSearch, "ref_long.txt"), LP);

            Save();
            if (ConsoleExport)
            {
                Console.WriteLine("Finished!");
                Console.ReadKey();
            }
            A = null;
            SL = null;
        }

        private static Aligner RunOptimizations(Aligner A, Seq_List SL)
        {
            //Test Conditions - Use this to make sure that your settings are correct
            Aligner.s.MaxSequences = 15000;
            Aligner.s.AttemptJoin = true; Aligner.s.ExportLocation = false; Aligner.s.QualScoringSkip = true;
            string Name;
            int TestMinReads2Target = S.AlignSize - S.PieceSize - 3;
            for (byte JoinRevTypes = 0; JoinRevTypes <= 3; JoinRevTypes++)
            {
                A = new Aligner(SL, S.PieceSize, S.PathFastQToSearch, S); //9 is ideal for most things
                A.Preferences = new Library_Preference();
                Aligner.s.JoinRevType = JoinRevTypes;
                A.AlignFolder(S.PathFastQToSearch, TestMinReads2Target, 1, true);
                Name = TestMinReads2Target + "-" + JoinRevTypes + ".txt";
                A.ExportResults(S.PathFastQToSearch + "res_" + Name, Name);
            }

            return A;
        }

        private static void SetupSchemas()
        {
            bool Skip = true;
            if (Skip) return;
            //string plasmid5 = "GGAAAGGACGAAACACCG"; string plasmid3 = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGC";
            //SL.Schema.AddStatic("GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGG");
            //SL.Schema.AddRedundantBlock(); SL.Schema.AddStatic(plasmid5);
            //SL.Schema.AddUnique(); SL.Schema.AddStatic(plasmid3); SL.Schema.AddRedundantBlock();
            {
                //LP.AddB("21", "1.txt");
                //LP.AddB("22", "2.txt");
                //LP.AddB("23", "3.txt");
                //LP.AddB("24", "4.txt");
                //LP.AddB("29", "1.txt");
                //LP.AddB("30", "2.txt");
                //LP.AddB("31", "3.txt");
                //LP.AddB("32", "4.txt");
            }
            {
                //LP.AddB("3A1", "MS760_Balce.txt");
                //LP.AddB("4A1", "MS000_McAllaster.txt");
                //LP.AddB("2A1", "MS764_Fpool.txt");
                //LP.AddB("1H1", "MS764_Fpool.txt");
                //LP.AddB("1A1", "MS765_Golgi.txt");

                //LP.AddB("E02", "ASX_1.txt");
                //LP.AddB("F02", "ASX_2.txt");
                //LP.AddB("G02", "ASX_3.txt");
                //LP.AddB("H02", "ASX_4.txt");

                //LP.AddB("E02", "ASX_Full.txt");
                //LP.AddB("F02", "ASX_Full.txt");
                //LP.AddB("G02", "ASX_Full.txt");
                //LP.AddB("H02", "ASX_Full.txt");
            }
            {
                /*
                //This is the style to use for an actually pooled crispr screen
                Library_Preference LP = new Library_Preference();
                LP.AddB("Miller-Plate09-A05-", "20171025Nisanth.txt"); //1
                LP.AddB("Miller-Plate10-A05-", "20171025Nisanth.txt"); //2
                LP.AddB("Miller-Plate11-A05-", "20171025Nisanth.txt"); //3
                LP.AddB("Miller-Plate10-H01-", "20171025Nisanth.txt"); //4

                Seq_List SL = new Seq_List();
                SL.Schema.AddUnique();
                */

                /* - - 11/2/2017
                LP.Add("MS841_McAllaster.txt", "Miller-Plate20-E02-"); //0
                LP.Add("MS760_Balce_gRNAs.txt", "Miller-Plate01-Hxx-"); //0
                LP.Add("McAllaster_All_Oligos.txt", "Miller-Plate01-Hyy-"); //1
                LP.Add("MS764_Final_Pool.txt", "Miller-Plate01-H12-"); //2
                LP.Add("MS765 Final_Golgi_gRNAs.txt", "Miller-Plate02-H12-"); //3
                //--"PCR v01"
                LP.Add("MS760_Balce_gRNAs.txt", "Miller-Plate01-H09-"); //0
                LP.Add("McAllaster_All_Oligos.txt", "Miller-Plate01-H10-"); //1
                LP.Add("MS764_Final_Pool.txt", "Miller-Plate01-H11-"); //2
                LP.Add("MS765 Final_Golgi_gRNAs.txt", "Miller-Plate01-H12-"); //3
                //--"Library v01"
                LP.Add("MS760_Balce_gRNAs.txt", "Miller-Plate01-H12-"); //0
                LP.Add("McAllaster_All_Oligos.txt", "Miller-Plate02-H12-"); //1
                LP.Add("MS764_Final_Pool.txt", "Miller-Plate03-H12-"); //2
                LP.Add("MS765 Final_Golgi_gRNAs.txt", "Miller-Plate04-H12-"); //3
                */
            }
        }

        private static int FinalAlignSize; //Align size is requested but the file size may be adjusted
        public static int MinReads2Target
        {
            get
            {
                int Min = 1 + FinalAlignSize - (S.PieceSize * (S.Mismatches_Allowed + 1));
                Min = Math.Max(1, Min);
                Min = Math.Min(FinalAlignSize - 1, Min);
                return Min;
            }
        }

        public static Store S; //This holds the parameters
        private static string SavePath = Path.Combine(Path.GetTempPath(), "Aligner_Settings.xml");

        public static void Load()
        {
            Load(SavePath);
        }

        public static Store Load(string LoadPath)
        {
            SavePath = LoadPath;
            Store s = Store.Load(SavePath);
            return s;
        }

        public static void Save()
        {
            try
            {
                string LocalSave = Path.Combine(S.PathReferenceList, DateTime.Now.ToString("yyyyMMdd_hhmm") + "_Settings.xml");
                Save(LocalSave);
            }
            catch
            {
                Save(SavePath);
            }
        }

        public static void Save(string SavePathThis)
        {
            if (S == null) S = new Store();
            Store s = S;

            s.Save(SavePathThis);
            Console.WriteLine("Settings saved..");
            Console.WriteLine(SavePathThis);
        }

        public static void Message(string Message, bool WaitForKey, bool DebugBreak, System.ComponentModel.BackgroundWorker BW = null)
        {
            if (BW != null) BW.ReportProgress(0, Message);
            if (_ConsoleExport)
            {
                Console.WriteLine(Message);
                if (WaitForKey) Console.ReadKey();
            }
            if (DebugBreak) System.Diagnostics.Debugger.Break();
        }
    }


    public class Store
    {
        public string PathReferenceList;
        public string PathFastQToSearch;

        public int PieceSize;
        public int AlignSize;

        public bool AttemptJoin;
        public byte JoinRevType;
        public bool FailedJoin_UseBestRead; //If the Join Fails, just use the read with the longest good quality score (TRUE), if false, will keep both reads, and append "NNN" in the middle

        public bool QualScoringSkip;

        public int Mismatches_Allowed;
        public int MaxSitesPerAmplicon;
        public bool ReverseComp_IncomingRefSequence;

        public int MinutesMaxPerIndex;
        public int MinSequencesForParallelAlignment;
        public int MaxSequences;
        public int Join_RegionSearch;
        public int Join_RegionSearchStep;
        public int Join_RegionSearchStart;
        public int Join_AttemptMatch_Longest;
        public int Join_AttemptMatch_Shortest;
        public int Join_AttemptMatch_Step;

        public bool ExportLocation;
        public bool Demux;

        public int NumberOfFastASeq;
        public int RefTopSeq_NumSeq;

        [System.ComponentModel.Description("Testing this out to get the right metadata")]
        public bool SplitLargeFiles;

        public static Store GetDefaults()
        {
            //System.Diagnostics.Debugger.Break(); //Change below as needed
            Store iS = new Store();

            iS.PathReferenceList = @"R:\dB\Sequences\gRNA Pools\";
            iS.PathFastQToSearch = @"c:\temp\ngs\ngs\";

            iS.SplitLargeFiles = false;
            iS.ReverseComp_IncomingRefSequence = false;

            // Important 8/16/2018 - Can't run in parallel if demuxing . . for now we have to only do direct
            iS.MinSequencesForParallelAlignment = 15000; //Keep this relatively high 15000, or make it way high to turn off parallelization
            iS.Join_RegionSearch = 105; //60 defualt
            iS.Join_RegionSearchStep = 4; //4 default
            iS.Join_RegionSearchStart = 4;
            iS.Join_AttemptMatch_Longest = 15; //15 default
            iS.Join_AttemptMatch_Shortest = 12; //12 default
            iS.Join_AttemptMatch_Step = 3;   //3 default
            iS.FailedJoin_UseBestRead = true; //Will just use one of the reads, not both if set to TRUE . . FALSE will use both with NNN stuck in between
            iS.NumberOfFastASeq = 250;       //Useful if you want to align the outputted sequences
            iS.RefTopSeq_NumSeq = 25;
            iS.MaxSitesPerAmplicon = 3; //For gRNA Libraries, this is set to 1 (there is only 1 gRNA per sequence), but for amplicons where there is more than one mutation site per amplicon, list the # of sites per amplicon

            iS.PieceSize = 20;  //Usually 7 With a 20 mer, you can use a piece size of 9 with 11 matches to get a perfect sequence.
                                //If you are allowing mismatches, then do something like 6 x 12
            iS.AlignSize = -1; // -1 //Usually 20, Asxl1 36
            iS.Mismatches_Allowed = 0; //Usually 3
            iS.MinutesMaxPerIndex = 7;

            //NGS CRISPR Library
            iS.Demux = false; //Can't do parallel alignment if demuxing7
            iS.ExportLocation = iS.Demux;
            iS.QualScoringSkip = true;
            iS.AttemptJoin = false;
            iS.JoinRevType = 2; //Usually 0 or 2
            iS.MaxSequences = -1; // 1000; //-1;

            return iS;
        }

        public static string DefaultSavePath { get => Path.Combine(Path.GetTempPath(), "LA_StoreSettings.xml"); }

        public static Store Load()
        {
            return Load(DefaultSavePath);
        }

        public static Store Load(string savePath)
        {
            FileInfo FI = new FileInfo(savePath);
            if (!FI.Exists) return Store.GetDefaults();
            using (var stream = System.IO.File.OpenRead(savePath))
            {
                var serializer = new XmlSerializer(typeof(Store));
                return (Store)serializer.Deserialize(stream);
            }
        }

        public void Save()
        {
            Save(DefaultSavePath);
        }

        public void Save(string savePath)
        {
            using (var writer = new System.IO.StreamWriter(savePath))
            {
                var serializer = new XmlSerializer(this.GetType());
                serializer.Serialize(writer, this);
                writer.Flush();
            }
        }

        public Store_ParamsList GetList()
        {
            Store_ParamsList List = new Store_ParamsList();
            foreach (System.Reflection.FieldInfo FI in this.GetType().GetFields())
            {
                if (!FI.Name.StartsWith("_"))
                {
                    KVP_PI KI = new KVP_PI(FI.Name, FI.GetValue(this).ToString());
                    List.Add(KI);
                }
            }
            return List;
        }

        public void UpdateFrom(IEnumerable<KVP_PI> dataSource)
        {
            Store_ParamsList NA_List = new Store_ParamsList(dataSource);
            Type t = this.GetType(); System.Reflection.FieldInfo FI;
            foreach (KVP_PI kI in NA_List)
            {
                FI = t.GetField(kI.Key);
                switch (FI.FieldType.Name)
                {
                    case "Byte":
                        FI.SetValue(this, byte.Parse(kI.Value)); break;
                    case "Int32":
                        FI.SetValue(this, int.Parse(kI.Value)); break;
                    case "Double":
                        FI.SetValue(this, double.Parse(kI.Value)); break;
                    case "Boolean":
                        FI.SetValue(this, bool.Parse(kI.Value)); break;
                    default:
                        FI.SetValue(this, kI.Value); break;
                }
            }
        }
    }

    public class KVP_PI
    {
        public string Key { get; }
        public string Value { get; set; }

        public KVP_PI(string k, string v)
        {
            Key = k; Value = v;
        }
    }

    public class Store_ParamsList : IEnumerable<KVP_PI>
    {
        List<KVP_PI> List;

        public Store_ParamsList()
        {
            List = new List<KVP_PI>();
        }

        public Store_ParamsList(IEnumerable<KVP_PI> dataSource)
        {
            List = dataSource.ToList();
        }

        public void Add(KVP_PI Item)
        {
            List.Add(Item);
        }

        public int Count => List.Count;

        public Store GetParams()
        {
            Store Params = new Store(); Type t = Params.GetType();
            System.Reflection.FieldInfo FI;
            foreach (KVP_PI kI in List)
            {
                FI = t.GetField(kI.Key);
                FI.SetValue(Params, kI.Value);
            }
            return Params;
        }

        public IEnumerator<KVP_PI> GetEnumerator()
        {
            return List.GetEnumerator();
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return List.GetEnumerator();
        }
    }

    public class Library_Preference
    {
        public Dictionary<string, short> ReferenceName_rIndex_Hash;
        public Dictionary<string, short> IndexName_rIndex_Hash;
        public List<string> ReferenceName;
        public List<string> IndexName_StartsWith;
        public List<string> UnMatched { get; internal set; }

        public Library_Preference()
        {
            ReferenceName_rIndex_Hash = new Dictionary<string, short>();
            IndexName_rIndex_Hash = new Dictionary<string, short>();
            ReferenceName = new List<string>();
            IndexName_StartsWith = new List<string>();
            UnMatched = new List<string>();
        }

        public void Add(string referenceName, string indexName)
        {
            short idx;
            idx = (short)ReferenceName.Count;
            ReferenceName.Add(referenceName);
            IndexName_StartsWith.Add(indexName);
            ReferenceName_rIndex_Hash.Add(referenceName, idx);
            IndexName_rIndex_Hash.Add(indexName, idx);
        }

        public void AddB(string NGSindexName, string LibraryreferenceName)
        {
            short idx;
            if (ReferenceName_rIndex_Hash.ContainsKey(LibraryreferenceName.ToUpper()))
            {
                idx = ReferenceName_rIndex_Hash[LibraryreferenceName.ToUpper()];
            }
            else
            {
                idx = (short)ReferenceName.Count;
                ReferenceName.Add(LibraryreferenceName.ToUpper());
                ReferenceName_rIndex_Hash.Add(LibraryreferenceName.ToUpper(), idx);
            }
            IndexName_StartsWith.Add(NGSindexName);
            IndexName_rIndex_Hash.Add(NGSindexName, idx);
        }

        public short Get_orAdd_ReferenceIndex(string FileName)
        {
            if (!ReferenceName_rIndex_Hash.ContainsKey(FileName.ToUpper()))
            {
                Add(FileName.ToUpper(), Path.GetFileNameWithoutExtension(FileName));
            }
            return ReferenceName_rIndex_Hash[FileName.ToUpper()];
        }

        public int Count { get => ReferenceName.Count; }

        public short GetIndexFrom_IndexFI(FileInfo FI)
        {
            if (Count <= 1)
            {
                return -1;
            }
            string Name = FI.Name.ToUpper();
            int length = IndexName_StartsWith[0].Length;
            string tName = Name.Substring(0, Math.Min(Name.Length, length));
            if (IndexName_rIndex_Hash.ContainsKey(tName))
            {
                return IndexName_rIndex_Hash[tName];
            }
            foreach (KeyValuePair<string, short> KVP in IndexName_rIndex_Hash)
            {
                if (Name.Contains(KVP.Key))
                {
                    return KVP.Value;
                }
            }
            return -1;
        }

        public string GetIndexShortName(FileInfo FI)
        {
            string Check = FI.Name.ToUpper();
            /*
            foreach (string ShortName in IndexName_StartsWith)
            {
                if (Check.Contains(ShortName.ToUpper()))
                {
                    return ShortName;
                }
            }
            */
            string t = Check.Replace("MILLER-PLATE", "");
            t = t.Replace("_001.FASTQ", "");
            t = t.Replace(".FASTQ", "");
            t = t.Replace("GEIC_PLATE", "");
            t = t.Replace("GEIC-PLATE", "");
            t = t.Replace("BUCHSER-PLATE", "");
            t = t.Replace("BUCHSER_PLATE", "");
            t = t.Replace(".ASSEMBLED", "");
            UnMatched.Add(t);
            return t;
        }
    }

    public class SourceStats
    {
        public int RefIndex = 0;
        public int Sequences = 0;
        public int Sequences_Unique = 0;
        public int Sequences_Failed = 0;
        public double Sequences_Mapped = 0;

        public SourceStats(int reft, int seq, int seq_uniq, int seq_failed)
        {
            RefIndex = reft; Sequences = seq; Sequences_Unique = seq_uniq; Sequences_Failed = seq_failed;
        }
    }

    public class Aligner
    {
        public Library_Preference Preferences;
        private Ref_Index _RefIndex;
        private List<string> _Results;
        public AlignedResultsManager AResultsManager;
        private StreamWriter _SeqTable_Exporter;
        public Seq_List RefList { get; private set; }
        public bool ExportSeqTable = false;
        public Dictionary<string, SourceStats> SourceStats { get; internal set; }
        private Dictionary<string, StreamWriter> DemuxWriters;
        private Dictionary<string, int> Track;
        private int FailTrack_Count;
        private int FailTrack_Total;

        public static Store s;

        public static int MinSequencesToAlign = 1; //Set this to 0 if you want to analyze files that have small # of sequences
        public static int Qual_MaxQScore = 26;
        public static int Qual_MinQScore = 10;
        public static int Qual_MinLength = 80;

        public static bool QualScoringParallel = true;

        public static bool Demux_Search_Main_RevComp = true;
        public static double Demux_Search_Start_Fract = 0.79;
        public static double Demux_Search_Stop_Fract = 0.99;
        public static bool Demux_Search_OtherSide = true;//Will also search the oppopsite on the other side

        public static StreamWriter TopFastASeq_SW;

        public Aligner(Seq_List SeqList, int BitSize, string FolderPathFastQ, Store store)
        {
            Aligner.s = store; //This isn't the best way, but at least we have one set of parameters now

            Ref_Index tRI = new Ref_Index(SeqList, BitSize);
            init(tRI, SeqList, FolderPathFastQ);
        }

        public Aligner(Ref_Index Ref, Seq_List SeqList, string FolderPathFastQ)
        {
            init(Ref, SeqList, FolderPathFastQ);
        }

        private void init(Ref_Index Ref, Seq_List SeqList, string FolderPathFastQ)
        {
            Preferences = new Library_Preference();
            _RefIndex = Ref;
            RefList = SeqList;
            _Results = new List<string>();
            AResultsManager = new AlignedResultsManager(FolderPathFastQ, SeqList);
            SourceStats = new Dictionary<string, SourceStats>();
        }

        private StreamWriter SWGeneral;
        private StreamWriter SWLocations;
        private long SWLocation_Seq_Index;

        public static bool SkipRanking = false;
        private static int NextNum = 0;

        public void DemuxFolder(string FolderPath, string DestPath, int MisMatchesAllowed)
        {
            //Mismatches not implemented
            DirectoryInfo DI = new DirectoryInfo(FolderPath);
            List<FileInfo> Files = GetSequencingFiles(DI);
            if (s.ExportLocation) SWLocations = new StreamWriter(Path.Combine(DestPath, "Locations.txt"));
            if (Files.Count > 0)
            {
                //First setup all the StreamWriters - - - - - - - - 
                DemuxWriters = new Dictionary<string, StreamWriter>(RefList.List.Count);
                Track = new Dictionary<string, int>(RefList.List.Count);
                foreach (Seq seq in RefList.List)
                {
                    DemuxWriters.Add(seq.Name, new StreamWriter(Path.Combine(DestPath, seq.Name + "_" + seq.OriginalSeq + ".fastq")));
                    Track.Add(seq.Name, 0);
                }

                //Now do the Demuxing - - - - - - - - - - - - - - - - -
                FailTrack_Count = 0; FailTrack_Total = 0; int idx = 0;
                foreach (FileInfo FI in Files) //Convert to Parallel
                {
                    DemuxFile(FI.FullName, MisMatchesAllowed);
                    if (idx++ % 50 == 0)
                        Console.WriteLine((100 * idx / Files.Count).ToString("00.0") + "% Done, " + (100 * FailTrack_Count / FailTrack_Total).ToString("00.0") + "% Failed");
                }

                //Now Close all the streamwriters - - - - - - - - - - -
                foreach (KeyValuePair<string, StreamWriter> KVP in DemuxWriters) KVP.Value.Close();
            }
            else
            {
                Console.WriteLine("No Sequence Files Found");
            }
            if (s.ExportLocation) SWLocations.Close();
            StringBuilder sb = new StringBuilder();
            foreach (KeyValuePair<string, int> KVP in Track) sb.Append(KVP.Key + "\t" + KVP.Value + "\r\n");
            File.WriteAllText("e:\\temp\\NGS\\demux_track.txt", sb.ToString());
        }

        public void OnlyJoin(string FolderPath)
        {
            DirectoryInfo DI = new DirectoryInfo(FolderPath);
            List<FileInfo> Files = GetSequencingFiles(DI);
            foreach (FileInfo FI in Files)
            {
                JoinFiles(FI.FullName);
            }
        }

        public void AlignFolder(string FolderPath, int MinReadsToTarget, int MinutesMax_PerFile, bool AllowMultipleMatches, System.ComponentModel.BackgroundWorker BW = null)
        {
            if (ExportSeqTable) _SeqTable_Exporter = new StreamWriter(FolderPath + "_SeqTableExport.txt");
            SWGeneral = new StreamWriter(Path.Combine(FolderPath, "Report" + (Aligner.NextNum++) + ".txt"));
            SWGeneral.WriteLine(DateTime.Now);
            SWGeneral.WriteLine(Aligner.s.QualScoringSkip.ToString() + " QualSkip, " + Aligner.s.JoinRevType.ToString() + " RevStandrd, " + Aligner.s.AttemptJoin.ToString() + " AttemptJoin");
            DirectoryInfo DI = new DirectoryInfo(FolderPath);
            short RefIndex = -1;
            string IndexShortName = "";
            List<FileInfo> Files = GetSequencingFiles(DI);
            if (Files.Count < 1)
            {
                Message("Found No FastQ files, or no Joinable files.", BW);
                return;
            }
            else
            {
                Message("Found " + Files.Count + " FastQ (pairs).", BW);
            }
            if (BW != null) if (BW.CancellationPending) return;
            foreach (FileInfo FI in Files)
            {
                RefIndex = Preferences.GetIndexFrom_IndexFI(FI);
                IndexShortName = Preferences.GetIndexShortName(FI);
                AlignFile(FI.FullName, false, RefIndex, IndexShortName, MinReadsToTarget, MinutesMax_PerFile, AllowMultipleMatches, s.NumberOfFastASeq > 0, BW);
                if (BW != null) if (BW.CancellationPending) break;
            }
            if (ExportSeqTable) _SeqTable_Exporter.Close();
            if (SWGeneral != null) SWGeneral.Close();
            if (FastaCompileSW != null) FastaCompileSW.Close(); FastaCompileSW = null;
            if (TopFastASeq_SW != null) TopFastASeq_SW.Close(); TopFastASeq_SW = null;
        }

        public static void Message(string Message, System.ComponentModel.BackgroundWorker BW = null)
        {
            if (BW == null)
                Console.WriteLine(Message);
            else
                BW.ReportProgress(0, Message);
        }

        private static List<FileInfo> GetSequencingFiles(DirectoryInfo DI)
        {
            //Searches several combinations in subfolders
            //AttemptJoin = false;
            List<string> SearchString = s.AttemptJoin ? 
                new List<string>(3) { "*_R1_0*.fastq*", "*_R1__*.fastq*", "*_R1.fastq*" } :
                new List<string>(3) { "*.fastq", "*.fq*", "*.fa" };
            List<FileInfo> Files = new List<FileInfo>();
            for (int i = 0; i < SearchString.Count; i++)
            {
                Files.AddRange(DI.GetFiles(SearchString[i]).ToList());
                foreach (DirectoryInfo DIs in DI.GetDirectories()) Files.AddRange(DIs.GetFiles(SearchString[i]));
                if (Files.Count > 1) break;
            }
            Files.RemoveAll(x => x.Name.StartsWith("."));
            return Files;
        }

        private void DemuxFile(string FilePath, int MisMatchesAllowed)
        {
            using (StreamReader SR = new StreamReader(FilePath))
            {
                string Seq; string Header; string Quality; string t;
                while (!SR.EndOfStream)
                {
                    t = SR.ReadLine();
                    if (t.StartsWith("@"))
                    {
                        Header = t;
                        Seq = SR.ReadLine();
                        if (SR.Peek() == 43)
                        {
                            SR.ReadLine();
                            Quality = SR.ReadLine();
                        }
                        else
                        {
                            Quality = "";
                        }
                        DemuxFileSeq(Header, Seq, Quality);
                    }
                }
                SR.Close();
            }
        }

        public string PrePostSeq(string Seq, string SeqBit, int index, int Amount, bool Pre)
        {
            int start;
            int end;
            if (Pre)
            {
                start = index - Amount;
                start = start < 0 ? 0 : start;
                end = index;
            }
            else
            {
                start = index + SeqBit.Length;
                end = start + Amount;
                end = end > Seq.Length ? Seq.Length : end;
            }
            return Seq.Substring(start, end - start);
        }

        private void DemuxFileSeq(string Header, string seq, string Quality)
        {
            SWLocation_Seq_Index++;
            FailTrack_Total++;
            List<Seq> Possible;
            string seqbit;
            string pre_rev = "";
            int start = (int)(seq.Length * (!Demux_Search_Main_RevComp ? Demux_Search_Start_Fract : (1 - Demux_Search_Stop_Fract))); // (int)(seq.Length * Demux_Search_Start_Fract);
            int end = (int)(seq.Length * (!Demux_Search_Main_RevComp ? Demux_Search_Stop_Fract : (1 - Demux_Search_Start_Fract))); // (int)(seq.Length * Demux_Search_Stop_Fract);
                                                                                                                                   //We also know for the usual ones that right after the barcode there is a N ATCTCGTAT
                                                                                                                                   //Look for the regular version in the beginning, but we need to invert the whole sequence
                                                                                                                                   //Look for the revcomp version at the end, and keep these like they are
            if (Demux_Search_Main_RevComp)
            {
                pre_rev = seq;
                seq = SeqHelper.ReverseComplement(seq);
            }
            for (int first = 0; (Demux_Search_OtherSide ? first < 2 : first < 1); first++)
            {
                if (first == 1)
                {
                    if (pre_rev == "")
                    {
                        pre_rev = seq;
                        seq = SeqHelper.ReverseComplement(seq);
                    }
                    else
                    {
                        string t = seq;
                        seq = pre_rev;
                        pre_rev = t;
                    }
                }
                for (int i = start; i < Math.Min(end, seq.Length - _RefIndex.BitSize); i++)
                {
                    seqbit = seq.Substring(i, _RefIndex.BitSize);
                    if (_RefIndex.MainHash.ContainsKey(seqbit))
                    {
                        Possible = _RefIndex.MainHash[seqbit];
                        Track[Possible[0].Name]++;
                        if (s.ExportLocation)
                        {
                            SWLocations.WriteLine(SWLocation_Seq_Index.ToString("X") + "\t" + seq.Length + "\t" + Possible[0].Name + "\t" +
                                i + "\t" + PrePostSeq(seq, seqbit, i, 10, true) + "\t" + PrePostSeq(seq, seqbit, i, 10, false));
                        }
                        else
                        {
                            //DemuxWriters[Possible[0].Name].WriteLine(Possible[0].Name + "_" + Possible[0].OriginalSeq);
                            DemuxWriters[Possible[0].Name].WriteLine(Header);
                            DemuxWriters[Possible[0].Name].WriteLine(Demux_Search_Main_RevComp ? pre_rev : seq);
                            DemuxWriters[Possible[0].Name].WriteLine("+");
                            DemuxWriters[Possible[0].Name].WriteLine(Quality);
                            return;
                        }
                    }
                }
            }
            //System.Diagnostics.Debug.Print(sB.ToString());
            FailTrack_Count++;
        }

        //start = (int)(seq.Length * (first == 0 ? Demux_Search_Start_Fract:(1- Demux_Search_Start_Fract)));
        //end = (int)(seq.Length* (first == 0 ? Demux_Search_Stop_Fract : (1 - Demux_Search_Stop_Fract)));

        private static StreamWriter FastaCompileSW = null;

        public void Export_ShortFASTA(List<string> Sequences, int Number, string FileName)
        {
            string FileNameNoExt = Path.GetFileNameWithoutExtension(FileName);
            if (FastaCompileSW == null)
            {
                string NewFileName = FileNameNoExt + ".FASTA";
                string Folder = Path.GetDirectoryName(FileName);
                string NewFolder = Path.Combine(Folder, "FastA");
                DirectoryInfo DI = new DirectoryInfo(NewFolder);
                if (!DI.Exists) DI.Create();
                FastaCompileSW = new StreamWriter(Path.Combine(NewFolder, NewFileName));
                for (int i = 0; i < RefList.List.Count; i++)
                {
                    FastaCompileSW.WriteLine(">Ref " + RefList.List[i].Name);
                    FastaCompileSW.WriteLine(RefList.List[i].OriginalSeq);
                }
            }
            for (int i = 0; i < Math.Min(Number, Sequences.Count); i++)
            {
                FastaCompileSW.WriteLine(">" + FileNameNoExt + " " + i.ToString("00000"));
                FastaCompileSW.WriteLine(Sequences[i]);
            }
        }

        public void JoinFiles(string FilePath)
        {
            int TotalCount; double JoinSuccessFraction;
            Dictionary<string, int> HashSequences_Ready;
            AlignFile_Pre(FilePath, false, null, true, out TotalCount, out HashSequences_Ready, out JoinSuccessFraction);
        }

        public void AlignFile(string FilePath, bool ReverseCompliment, short ReferenceIndex, string IndexShortName, int MinReadsToTarget, int MinutesMax, bool AllowMultipleMatches, bool ExportShortFastA, System.ComponentModel.BackgroundWorker BW = null)
        {
            //Doing this separately so that the temporary lists can be discarded
            int TotalCount; int Counter = 0;  int Failed = 0; double JoinSuccessFraction;
            Dictionary<string, int> HashSequences_Ready; FileInfo FI = new FileInfo(FilePath);
            Message("  Loading " + IndexShortName + ", " + (FI.Length/1024/1024).ToString("0.0") + " Mb.", BW); if (BW != null) if (BW.CancellationPending) return;
            AlignFile_Pre(FilePath, ReverseCompliment, ExportSeqTable ? _SeqTable_Exporter : null, out TotalCount, out HashSequences_Ready, out JoinSuccessFraction);

            Message("  Prep Aligning . . " + HashSequences_Ready.Count.ToString("0,000") + " Unique Sequences.", BW); if (BW!=null) if (BW.CancellationPending) return;
            string FileName = FI.Name;
            List<string> Sequences = HashSequences_Ready.Keys.ToList();

            if (SourceStats.Count == 0) SourceStats.Add("", new SourceStats(0, 0, 0, 0));
            SourceStats.Add(IndexShortName, new SourceStats(ReferenceIndex, TotalCount, Sequences.Count, 0));

            if (ExportShortFastA) Export_ShortFASTA(Sequences, s.NumberOfFastASeq, FilePath);
            Export_TopSequences(HashSequences_Ready, FilePath, IndexShortName, TotalCount);
            DateTime tStart = DateTime.Now; bool GaveUpdate = false;
            void CheckGiveUpdate(float UnMapped)
            {
                if (DateTime.Now > tStart.AddSeconds(20))
                {
                    GaveUpdate = true;
                    Message("  Getting " + UnMapped.ToString("0.0%") + " UnMapped Reads so Far . . ", BW);
                }
            }
            try
            {
                if (HashSequences_Ready.Count > MinSequencesToAlign)
                {
                    bool Direct = HashSequences_Ready.Count < s.MinSequencesForParallelAlignment; //Not worth setting up the parallel unless it is really big
                    Message("  Aligning with " + (Direct ? "Single Thread" : "Multi Threads"), BW);
                    if (Direct)
                    {
                        foreach (string sequence in Sequences)
                        {
                            Failed += AlignSeq(sequence, IndexShortName, HashSequences_Ready[sequence], ReferenceIndex, MinReadsToTarget, AllowMultipleMatches) ? 0 : 1;
                            if (BW != null) if (BW.CancellationPending) break;
                            Counter++;
                            if (!GaveUpdate) CheckGiveUpdate((float)Failed/Counter);
                        }
                    }
                    else
                    {
                        Parallel.ForEach(Sequences, (sequence, state) =>
                        {
                            if (!AlignSeq(sequence, IndexShortName, HashSequences_Ready[sequence], ReferenceIndex, MinReadsToTarget, AllowMultipleMatches))
                                Interlocked.Add(ref Failed, 1);
                            if (DateTime.Now > tStart.AddMinutes(MinutesMax)) state.Break(); //Limit this to three minutes
                            if (BW != null) if (BW.CancellationPending) state.Break();
                            Interlocked.Add(ref Counter, 1);
                            if (!GaveUpdate) CheckGiveUpdate((float)Failed / Counter);
                        });
                    }
                }
            }
            catch (Exception E)
            {
                Console.WriteLine("Aligner Error: " + E.Message);
                try
                {
                    Console.ReadKey();
                    Environment.Exit(-1);
                }
                catch { }
            }
            SourceStats[IndexShortName].Sequences_Failed = Failed;
            string Report = TotalCount.ToString() + " Reads, " + Sequences.Count.ToString() + " Unique, " + ((double)Failed / Sequences.Count).ToString("0.0%") + " UnMapped, " + (double.IsNaN(JoinSuccessFraction) ? "" : JoinSuccessFraction.ToString("0.0%") + " Joined");
            Message(Report, BW);
            SWGeneral.WriteLine(Report);

            HashSequences_Ready = null;
            Sequences = null;
        }

        private void Export_TopSequences(Dictionary<string, int> hashSequences_Ready, string filePath, string index, int TotalCount)
        {
            //First build the sorted list . . don't bother if there are less than 2 sequences observed
            SortedList<double, string> S = new SortedList<double, string>(); double Key;
            double Max = 0; double Threshold = 2;
            foreach (KeyValuePair<string, int> keyValuePair in hashSequences_Ready)
            {
                if (keyValuePair.Value >= Threshold)
                {
                    Key = -keyValuePair.Value;
                    while (S.ContainsKey(Key)) Key -= 0.0000001;
                    S.Add(Key, keyValuePair.Key);
                    if (Max < -Key)
                    {
                        Max = -Key; Threshold = Math.Max(2, Max / 10);
                    }
                }
            }
            // Now restrict the list to the Top X
            if (TopFastASeq_SW == null)
            {
                TopFastASeq_SW = new StreamWriter(Path.Combine(Path.GetDirectoryName(filePath), "ref_TopSeq.txt"));
                TopFastASeq_SW.WriteLine("Index\tTotal\tCount\tSequence");
            }
            for (int i = 0; i < Math.Min(s.RefTopSeq_NumSeq, S.Count); i++)
            {
                TopFastASeq_SW.WriteLine(index + "\t" + TotalCount + "\t" + (-S.Keys[i]) + "\t" + S.Values[i]);
            }
        }

        //So it doesn't break older versions
        private static void AlignFile_Pre(string FilePath, bool ReverseCompliment, StreamWriter SW, out int TotalSequenceCount, out Dictionary<string, int> HashSequences_Ready, out double JoinSuccessFraction)
        {
            AlignFile_Pre(FilePath, ReverseCompliment, SW, false, out TotalSequenceCount, out HashSequences_Ready, out JoinSuccessFraction);
        }

        private static void AlignFile_Pre(string FilePath, bool ReverseCompliment, StreamWriter SW, bool OnlyWriteOutNewJoinedReads, out int TotalSequenceCount, out Dictionary<string, int> HashSequences_Ready, out double JoinSuccessFraction)
        {
            Console.WriteLine("");
            Console.WriteLine("Reading " + FilePath);
            ConcurrentBag<string> ListSequences_Quality;
            //Dictionary<string, string> SeqNames;
            if (FilePath.Contains(".faa") || FilePath.Contains(".fasta"))
            {
                //FASTA Doesn't have Quality Scores
                ListSequences_Quality = AlignFile_Load_FASTA(FilePath, SW);
                JoinSuccessFraction = -1;
            }
            else
            {
                //FASTQ Needs to be loaded differently and has Quality Scores
                Tuple<List<List<string>>, double> tRet = AlignFile_Load_FASTQ_AttemptJoin(FilePath, s.AttemptJoin);
                List<List<string>> ListSequences_ReadIn = tRet.Item1;
                JoinSuccessFraction = tRet.Item2;
                if (OnlyWriteOutNewJoinedReads)
                {
                    OnlyWriteOutNewJoinedReadsNow(FilePath, JoinSuccessFraction, ListSequences_ReadIn);
                    TotalSequenceCount = ListSequences_ReadIn.Count;
                    HashSequences_Ready = new Dictionary<string, int>();
                    return;
                }
                ListSequences_Quality = new ConcurrentBag<string>();
                if (s.QualScoringSkip)
                {
                    foreach (List<string> seq_quality in ListSequences_ReadIn) ListSequences_Quality.Add(seq_quality[0]);
                }
                else
                {
                    Console.WriteLine("  Quality Window . . .  ");
                    if (QualScoringParallel)
                    {
                        Parallel.ForEach(ListSequences_ReadIn, (seq_quality) =>
                        {
                            string seq2 = GoodQualSequence_MinMax(seq_quality[0], seq_quality[1], Qual_MaxQScore, Qual_MinQScore, Qual_MinLength);
                            if (ReverseCompliment) seq2 = SeqHelper.ReverseComplement(seq2);
                            ListSequences_Quality.Add(seq2);
                        });
                    }
                    else
                    {
                        string seq2;
                        foreach (List<string> seq_quality in ListSequences_ReadIn)
                        {
                            seq2 = GoodQualSequence_MinMax(seq_quality[0], seq_quality[1], Qual_MaxQScore, Qual_MinQScore, Qual_MinLength);
                            if (ReverseCompliment) seq2 = SeqHelper.ReverseComplement(seq2);
                            ListSequences_Quality.Add(seq2);
                        }
                    }
                }
            }
            TotalSequenceCount = ListSequences_Quality.Count;

            Console.WriteLine("  Consolidating " + String.Format("{0:n0}", ListSequences_Quality.Count) + " Sequences . . . . ");
            HashSequences_Ready = new Dictionary<string, int>();
            foreach (string seq3 in ListSequences_Quality)
            {
                if (seq3 == null) continue;
                if (!HashSequences_Ready.ContainsKey(seq3))
                    HashSequences_Ready.Add(seq3, 0);
                HashSequences_Ready[seq3]++;
            }
        }

        private static void OnlyWriteOutNewJoinedReadsNow(string filePath, double joinSuccessFraction, List<List<string>> listSequences_ReadIn)
        {
            FileInfo FI = new FileInfo(filePath);
            DirectoryInfo DI = new DirectoryInfo(Path.Combine(FI.Directory.FullName, "JoinedRT" + Aligner.s.JoinRevType));
            if (!DI.Exists) DI.Create();
            string NewName = Path.Combine(DI.FullName, Path.GetFileNameWithoutExtension(FI.Name) + (joinSuccessFraction * 100).ToString(" 00") + Path.GetExtension(FI.Name));
            using (StreamWriter SW = new StreamWriter(NewName))
            {
                for (int i = 0; i < listSequences_ReadIn.Count; i++)
                {
                    SW.WriteLine("@" + i);
                    SW.WriteLine(listSequences_ReadIn[i][0]);
                    SW.WriteLine("+");
                    SW.WriteLine(listSequences_ReadIn[i][1]);
                }
                SW.Close();
            }
        }

        private static ConcurrentBag<string> AlignFile_Load_FASTA(string FilePath, StreamWriter SW)
        {
            ConcurrentBag<string> ListSequences_ReadIn;
            string t, h, more;
            FileInfo FI = new FileInfo(FilePath);
            StreamReader SR = new StreamReader(FilePath);
            ListSequences_ReadIn = new ConcurrentBag<string>();
            t = SR.ReadLine();
            while (!SR.EndOfStream)
            {
                if (t.StartsWith(">"))
                {
                    h = t.Substring(1);
                    t = SR.ReadLine(); more = "";
                    while (!t.StartsWith(">"))
                    {
                        more += t;
                        if (SR.EndOfStream) break;
                        t = SR.ReadLine();
                    }
                    ListSequences_ReadIn.Add(more);
                    if (Aligner.s.MaxSequences > 1)
                    {
                        if (ListSequences_ReadIn.Count > Aligner.s.MaxSequences) break;
                    }
                    if (SW != null) SW.WriteLine(FI.Name + "\t" + h + "\t" + more);
                    //t = SR.ReadLine();
                }
                else
                {

                }
            }
            SR.Close();
            return ListSequences_ReadIn;
        }

        private static Tuple<List<List<string>>, double> AlignFile_Load_FASTQ_AttemptJoin(string FilePath, bool AttemptJoin)
        {
            Tuple<List<List<string>>, double> FinalReturn;
            List<List<string>> ListSequences_ReadIn;
            ListSequences_ReadIn = new List<List<string>>();
            string seq;
            string quality;
            string OtherFile = AttemptJoin ? Find_Cognate_Reads_File(FilePath) : "";
            double JoinSuccess = 0; double JoinTotal = 0;
            Tuple<List<string>, bool> JoinReturn;
            if (OtherFile == "")
            {
                StreamReader SR = new StreamReader(FilePath);
                while (!SR.EndOfStream)
                {
                    SR.ReadLine(); seq = SR.ReadLine();
                    SR.ReadLine(); quality = SR.ReadLine();
                    ListSequences_ReadIn.Add(new List<string>(2) { seq, quality });
                    if (Aligner.s.MaxSequences > 1)
                    {
                        if (ListSequences_ReadIn.Count > Aligner.s.MaxSequences) break;
                    }
                }
                SR.Close();
            }
            else
            {
                StreamReader SR1 = new StreamReader(FilePath);
                StreamReader SR2 = new StreamReader(OtherFile);
                string R1, R2, H1, H2, Q1, Q2;
                JoinReadsMismatch = 0;
                while (!SR1.EndOfStream)
                {
                    H1 = SR1.ReadLine(); H2 = SR2.ReadLine();
                    switch (Aligner.s.JoinRevType)
                    {
                        default: //AKA 0
                                 //This is somehow backwards from correct
                            R1 = SR1.ReadLine(); R2 = SeqHelper.ReverseComplement(SR2.ReadLine());
                            SR1.ReadLine(); SR2.ReadLine();
                            Q1 = SR1.ReadLine(); Q2 = new String(SR2.ReadLine().ToCharArray().Reverse().ToArray());
                            break;
                        case 1:
                            R1 = SeqHelper.ReverseComplement(SR1.ReadLine()); R2 = SR2.ReadLine();
                            SR1.ReadLine(); SR2.ReadLine();
                            Q1 = new String(SR1.ReadLine().ToCharArray().Reverse().ToArray()); Q2 = SR2.ReadLine();
                            break;
                        case 2:
                            //This is usually the correct one
                            R2 = SR1.ReadLine(); R1 = SeqHelper.ReverseComplement(SR2.ReadLine());
                            SR1.ReadLine(); SR2.ReadLine();
                            Q2 = SR1.ReadLine(); Q1 = new String(SR2.ReadLine().ToCharArray().Reverse().ToArray());
                            break;
                        case 3:
                            //Fully Reversed
                            R2 = SeqHelper.ReverseComplement(SR1.ReadLine()); R1 = SR2.ReadLine();
                            SR1.ReadLine(); SR2.ReadLine();
                            Q2 = new String(SR1.ReadLine().ToCharArray().Reverse().ToArray()); Q1 = SR2.ReadLine();
                            break;
                    }
                    JoinReturn = Join_Reads_L(H1, H2, R1, R2, Q1, Q2); //ListSequences_ReadIn.Add(Join_Reads(H1, H2, R1, R2, Q1, Q2));
                    if (JoinReturn == null) {
                        //something went wrong, could be the end of one of the files, return here . . 
                        break;
                    }
                    ListSequences_ReadIn.Add(JoinReturn.Item1);
                    if (Aligner.s.MaxSequences > 1)
                    {
                        if (ListSequences_ReadIn.Count > Aligner.s.MaxSequences) break;
                    }
                    JoinTotal++;
                    if (JoinReturn.Item2) JoinSuccess++;
                }
                SR1.Close(); SR2.Close();
            }
            FinalReturn = new Tuple<List<List<string>>, double>(ListSequences_ReadIn, JoinSuccess / JoinTotal);
            return FinalReturn;
        }

        public static int JoinReadsMismatch;

        public static Tuple<List<string>, bool> Join_Reads_L(string H1, string H2, string R1, string R2, string Q1, string Q2)
        {
            if (s.Join_RegionSearchStart > R1.Length + 2) s.Join_RegionSearchStart = (R1.Length / 4);
            if (s.Join_RegionSearch > R1.Length + 2) s.Join_RegionSearch = (R1.Length / 2);
            List<string> Results = new List<string>(3) { "", "", H1 };
            Tuple<List<string>, bool> Final;
            if (!H1.Contains(" ") || !H2.Contains(" ")) return null;
            if (H1.Substring(0, H1.IndexOf(" ") - 1) != H2.Substring(0, H2.IndexOf(" ") - 1))
            {
                JoinReadsMismatch++;
                //System.Diagnostics.Debugger.Break();
                //This means that the headers were different, so it may be the wrong cognate file.
            }
            for (int start = s.Join_RegionSearchStart; start < s.Join_RegionSearch; start += s.Join_RegionSearchStep) //Now try to find the point of overlap in the two halves by starting at 4 and working forward
            {
                for (int length = s.Join_AttemptMatch_Longest; length >= s.Join_AttemptMatch_Shortest; length -= s.Join_AttemptMatch_Step)
                {
                    if (R2.Contains(R1.Substring(start, length)))
                    {
                        int l1 = R2.IndexOf(R1.Substring(start, length));
                        Results = Join_Reads_Quality_Based(R1, R2, Q1, Q2, start - l1);
                        Final = new Tuple<List<string>, bool>(Results, true);
                        return Final;
                    }
                }
            }
            //If we got here, that means we failed at joining
            if (s.FailedJoin_UseBestRead)
            {
                bool Use1 = ScoreQualityLength(Q1) > ScoreQualityLength(Q2);
                Results[0] = Use1 ? R1 : R2;
                Results[1] = Use1 ? Q1 : Q2;
            }
            else
            {
                Results[0] = R1 + "NNN" + R2;
                Results[1] = Q1 + "AAA" + Q2;
            }
            Final = new Tuple<List<string>, bool>(Results, false);
            return Final;
        }

        public static double ScoreQualityLength(string QScore)
        {
            double score = 0;
            char[] cArr = QScore.ToCharArray();
            foreach (char c in cArr)
            {
                score += Math.Pow(c, 1.2);
            }
            return score;
        }

        public static List<string> Join_Reads(string H1, string H2, string R1, string R2, string Q1, string Q2)
        {
            JoinReadsMismatch = 0;
            Tuple<List<string>, bool> Final = Join_Reads_L(H1, H2, R1, R2, Q1, Q2);
            return Final.Item1;
        }

        public static List<string> Join_Reads_Quality_Based(string R1, string R2, string Q1, string Q2, int StartDiff)
        {
            bool R2Before = StartDiff >= 0;
            StartDiff = Math.Abs(StartDiff);
            string a1S = R2Before ? R1 : R2;
            string a2S = R2Before ? R2 : R1;
            string a1Q = R2Before ? Q1 : Q2;
            string a2Q = R2Before ? Q2 : Q1;
            StringBuilder sSequ = new StringBuilder();
            StringBuilder sQual = new StringBuilder();

            try
            {
                //First make the Pre Sequence - - - - - - - - - - - - - - - 
                if (StartDiff > 0)
                {
                    sSequ.Append(a1S.Substring(0, StartDiff));
                    sQual.Append(a1Q.Substring(0, StartDiff));
                }

                //Now Make the Overlapping Part - - - -  - - - - - - - - - - -
                int Overlap = Math.Min(a1S.Length - StartDiff, a2S.Length);
                char[] c1S = new char[Overlap]; char[] c1Q = new char[Overlap]; char[] c2S = new char[Overlap]; char[] c2Q = new char[Overlap];
                a1S.CopyTo(StartDiff, c1S, 0, Overlap);
                a1Q.CopyTo(StartDiff, c1Q, 0, Overlap);
                a2S.CopyTo(0, c2S, 0, Overlap);
                a2Q.CopyTo(0, c2Q, 0, Overlap);
                bool Score;
                for (int i = 0; i < c1S.Length; i++)
                {
                    Score = ((c1S[i] == 'N') ? false : ((c2S[i] == 'N') ? true : c1Q[i] >= c2Q[i]));
                    //System.Diagnostics.Debug.WriteLine(((int)c1Q[i]).ToString("00") + " " + ((int)c2Q[i]).ToString("00") + " | " + c1Q[i] + " " + c2Q[i] + " | " + c1S[i] + " " + c2S[i] + " = " + (Score ? c1S[i] : c2S[i]));
                    sSequ.Append(Score ? c1S[i] : c2S[i]);
                    sQual.Append(Score ? c1Q[i] : c2Q[i]);
                }

                //Now do the tail  - - - - - - - - - - - - - - - - - - - - - - 
                int Leftover = Math.Max(a1S.Length - StartDiff, a2S.Length) - Overlap;
                if (Leftover > 0)
                {
                    bool Overlapped = (a1S.Length - StartDiff) > a2S.Length;
                    sSequ.Append(Overlapped ? a1S.Substring(Overlap, Leftover) : a2S.Substring(Overlap, Leftover));
                    sQual.Append(Overlapped ? a1Q.Substring(Overlap, Leftover) : a2Q.Substring(Overlap, Leftover));
                }
            }
            catch
            {
                //If we encounter an error, just return the first read
                sSequ = new StringBuilder(R1); sQual = new StringBuilder(Q1);
            }
            return new List<string>(2) { sSequ.ToString(), sQual.ToString() };
        }

        public static string Find_Cognate_Reads_File(string FileToCheck)
        {
            FileInfo FI = new FileInfo(FileToCheck);
            DirectoryInfo DI = FI.Directory;
            string OtherFileName = FI.Name.Replace("_R1_0", "_R2_0");
            FileInfo FI2 = new FileInfo(Path.Combine(DI.FullName, OtherFileName));
            if (FI2.FullName == FI.FullName || !FI2.Exists)
            {
                OtherFileName = FI.Name.Replace("_R1_", "_R2_");
                FI2 = new FileInfo(Path.Combine(DI.FullName, OtherFileName));
            }
            return (FI2.Exists && FI2.FullName != FI.FullName) ? FI2.FullName : "";
        }

        public static string GoodQualSequence_MinMax(string seq, string quality, int Max_QScore, int Min_QScore, int MinLength)
        {
            int point = seq.Length > 120 ? 91 : seq.Length / 2;
            string tSeq;
            if (quality == "") return seq; //This means that it is assumed to be perfect sequence, no need to do windowing
            tSeq = GoodQualSequence(seq, quality, Max_QScore, point);
            if (tSeq.Length < MinLength)
            {
                int MaxLength = 0;
                string BestSeq = tSeq;
                for (int MQScore = Max_QScore; MQScore >= (Min_QScore); MQScore -= 2)
                {
                    for (point = 14; point < seq.Length; point += (MinLength / 2))
                    {
                        tSeq = GoodQualSequence(seq, quality, MQScore, point);
                        if (tSeq.Length > MaxLength)
                        {
                            BestSeq = tSeq;
                            MaxLength = tSeq.Length;
                        }
                    }
                    tSeq = BestSeq;
                    if (tSeq.Length >= MinLength) break;
                }
            }
            return tSeq;
        }

        public static string GoodQualSequence(string seq, string quality, int Min_QScore, int point)
        {
            byte[] qual;
            qual = Encoding.ASCII.GetBytes(quality);
            int start, end;
            for (start = point; start > 2; start--)
                if (qual[start] < (Min_QScore + 33)) break;
            for (end = point; end < seq.Length - 10; end++)
                if (qual[end] < (Min_QScore + 33)) break;
            return seq.Substring(start + 1, end - start);
        }

        public bool AlignSeq(string seq, string Prefix, int occurences, short ReferenceIndex, int MinReadsToTarget, bool AllowMultipleMatches)
        {
            //AllowMultiple doesn't do anything in this version
            Dictionary<Seq, int> Possibilities = new Dictionary<Seq, int>();
            Dictionary<Seq, int> Alternate_Possibilities = new Dictionary<Seq, int>();
            Dictionary<string, byte> Used = new Dictionary<string, byte>();
            string seqbit; int Max = 0; int RepresentsAllCount = 0;
            for (int i = 0; i < seq.Length - _RefIndex.BitSize; i++)
            {
                seqbit = seq.Substring(i, _RefIndex.BitSize);
                if (_RefIndex.MainHash.ContainsKey(seqbit))
                {
                    List<Seq> PossibleMatches = _RefIndex.MainHash[seqbit];
                    if (_RefIndex.TotalSequence_Count > 2 && PossibleMatches.Count > _RefIndex.TotalSequence_Count * 0.9)
                        RepresentsAllCount++;
                    else if (_RefIndex.TotalSequence_Count < 2 || PossibleMatches.Count < _RefIndex.SmallestLibraryCount)
                    {
                        foreach (Seq s in PossibleMatches)
                        {
                            if (s.RepresentsAll)
                                RepresentsAllCount++;
                            else
                            {
                                if ((PossibleMatches.Count < 2) || (PossibleMatches.Count < (_RefIndex.SmallestLibraryCount * 0.1)))
                                {
                                    if (!Possibilities.ContainsKey(s))
                                        Possibilities.Add(s, 0);
                                    Possibilities[s]++;
                                    if (Possibilities[s] > Max) Max = Possibilities[s];
                                }
                                else
                                {
                                    if (!Alternate_Possibilities.ContainsKey(s))
                                        Alternate_Possibilities.Add(s, 0);
                                    Alternate_Possibilities[s]++;
                                }
                            }
                        }
                    }
                } //else it isn't in the hash, so ignore
            }
            AlignResultsClass ARC = new AlignResultsClass(Prefix, occurences, seq, seq.Length, RepresentsAllCount, Max);
            if (Alternate_Possibilities.Count > 0)
            {
                List<Seq> tKeys = Possibilities.Keys.ToList();
                foreach (Seq Key in tKeys)
                {
                    if (Alternate_Possibilities.ContainsKey(Key))
                        Possibilities[Key] += Alternate_Possibilities[Key];
                    if (Max < Possibilities[Key])
                        Max = Possibilities[Key];
                }
            }
            if (Max >= MinReadsToTarget)
            {
                //Sort the possible results from best to worst
                SortedList<double, Seq> sorted = new SortedList<double, Seq>(); double d; Random R = new Random();
                foreach (KeyValuePair<Seq, int> KVP in Possibilities)
                {
                    if (KVP.Value > (Max / 1.45))
                    {
                        d = -KVP.Value - (R.NextDouble() / 10);
                        while (sorted.ContainsKey(d))
                            d -= (R.NextDouble() / 100);
                        sorted.Add(d, KVP.Key);
                    }
                }
                //Keep only the requested number of "sites"
                List<Seq> ResultsList = sorted.Values.ToList();
                if (Aligner.s.MaxSitesPerAmplicon == 1)
                {
                    if (ResultsList.Count > 1) ResultsList = AlignToSmallSet(seq, ResultsList, ReferenceIndex);
                }
                else if (ResultsList.Count > Aligner.s.MaxSitesPerAmplicon)
                {
                    ResultsList.RemoveRange(Aligner.s.MaxSitesPerAmplicon, ResultsList.Count - Aligner.s.MaxSitesPerAmplicon);
                }
                if (s.ExportLocation) ARC.Location = Location(ResultsList[0], seq);
                for (int i = 0; i < ResultsList.Count; i++)
                {
                    ARC.TargetList.Add(ResultsList[i].Name);
                    ResultsList[i].AddResults(ReferenceIndex,
                        ((Math.Min(1, (double)Max / (MinReadsToTarget - 0)) * (double)occurences) / ResultsList.Count),
                        Prefix);
                }
                SourceStats[Prefix].Sequences_Mapped += occurences;
            }
            ProcessResults(ARC);
            return Max >= MinReadsToTarget;
        }

        public int Location(Seq SequenceToFind, string WithinSequence)
        {
            string SeqSearch;
            int len = SequenceToFind.OriginalSeq.Length;
            int stop = (len / 2) - 4;
            for (int start = 0; start < stop; start++)
            {
                SeqSearch = SequenceToFind.OriginalSeq.Substring(start, len -= 2);
                if (WithinSequence.Contains(SeqSearch))
                {
                    return WithinSequence.IndexOf(SeqSearch);
                }
            }
            return 0;
        }

        public void ProcessResults(AlignResultsClass aRC)
        {
            AResultsManager.ExportResults(aRC);
            //_Results.Add(aRC.ExportOld()); //Old Way
        }

        public List<Seq> AlignToSmallSet(string SeqToAlign, List<Seq> Reference, short ReferenceIndex)
        {
            return AlignToSmallSet(SeqToAlign, Reference, ReferenceIndex, Math.Min(_RefIndex.BitSize, 7));
        }

        public List<Seq> AlignToSmallSet(string SeqToAlign, List<Seq> Reference, short ReferenceIndex, int BitSize)
        {
            PriorityRemove(Reference); if (Reference.Count == 1) return Reference; //This looks if there are preferential sequences and ones that can be removed
            List<double> Track = new List<double>(Reference.Count);
            string tSeq;
            for (int i = 0; i < SeqToAlign.Length - BitSize; i++)
            {
                tSeq = SeqToAlign.Substring(i, BitSize);
                for (int j = 0; j < Reference.Count; j++)
                {
                    if (i == 0) Track.Add(0);
                    if (Reference[j].OriginalSeq.Contains(tSeq)) Track[j]++;
                }
            }
            SortedList<double, List<Seq>> sorted = new SortedList<double, List<Seq>>();
            double fraction;
            for (int j = 0; j < Reference.Count; j++)
            {
                fraction = -Track[j] / (Math.Min(SeqToAlign.Length, Reference[j].OriginalSeq.Length) - BitSize);
                if (!sorted.ContainsKey(fraction))
                    sorted.Add(fraction, new List<Seq>());
                sorted[fraction].Add(Reference[j]);
            }

            if (ReferenceIndex > -1)
            {
                //Now go through and get rid of things that don't match the library (if ReferenceIndex = -1, then it isn't known what we are looking at)
                for (int i = 0; i < sorted.Values.Count; i++)
                {
                    List<Seq> ToRemove = new List<Seq>();
                    int GoodIndexCount = 0;
                    foreach (Seq s in sorted.Values[i])
                        if (s.Reference_Index != ReferenceIndex)
                            ToRemove.Add(s);
                        else
                            GoodIndexCount++;
                    if (GoodIndexCount > 0)
                    {
                        foreach (Seq s in ToRemove)
                            sorted.Values[i].Remove(s);
                        return sorted.Values[i];
                    }
                    //This is a backup, if there is relatively good match for the second set, then we can see if a correct match is lurking here
                }
            }
            return sorted.Values[0];
        }

        private void PriorityRemove(List<Seq> reference)
        {
            HashSet<Seq> ToRemove = new HashSet<Seq>();
            foreach (Seq item in reference)
            {
                if (item.Name.Contains("|")) ToRemove.Add(item);
            }
            foreach (Seq item in ToRemove)
            {
                reference.Remove(item);
                if (reference.Count == 1) return;
            }
        }

        public void ExportResults(string FullName)
        {
            //Results are actually already expored, we just have to close off the Files
            AResultsManager.CloseFiles();
            //ExportResults(FullName, ""); //Old Way
        }

        public void ExportResults(string FullName, string Prefix)
        {
            StreamWriter SW = new StreamWriter(FullName);
            SW.WriteLine("Pre" + "\t" + "Prefix" + "\t" + "Occurences" + "\t" + "Seq" + "\t" + "Good Qual Length" + "\t" + "Common Matches" + "\t" + "Reads to Target" + "\t" + "Target" + "\t" + "Target" + "\t" + "Target" + "\t" + "Target" + "\t" + "Target");
            for (int i = 0; i < _Results.Count; i++)
            {
                SW.WriteLine(Prefix + "\t" + _Results[i]);
            }
            SW.Close();
        }

        internal void SplitLargeFiles(string pathMatch, int MegaBytesMax)
        {
            DirectoryInfo DI = new DirectoryInfo(pathMatch);
            foreach (FileInfo file in DI.GetFiles())
            {
                if (file.Length / 1024 / 1024 > MegaBytesMax)
                {
                    SplitLargeFile(file.FullName, MegaBytesMax);
                }
            }
        }

        private void SplitLargeFile(string fullName, int megaBytesMax)
        {
            string ext; string t; int i = 0;
            if (fullName.EndsWith(".SKIP")) return;
            FileInfo FI = new FileInfo(fullName);
            int Lines = 3000 * megaBytesMax;
            ext = Path.GetExtension(fullName).Trim().ToUpper();
            Dictionary<string, char> ExtToBreakPoint = new Dictionary<string, char>() { { ".FASTQ", '@' }, { ".FASTA", '>' } };
            using (StreamReader SR = new StreamReader(fullName))
            {
                while (!SR.EndOfStream)
                {
                    StreamWriter SW = new StreamWriter(Path.Combine(FI.Directory.FullName, Path.GetFileNameWithoutExtension(fullName) + "__" + i++ + "" + ext));
                    for (int j = 0; j < Lines; j++)
                    {
                        SW.WriteLine(SR.ReadLine());
                        if (SR.EndOfStream) break;
                    }
                    if (SR.EndOfStream) break;
                    while (SR.Peek() != ExtToBreakPoint[ext])
                    {
                        t = SR.ReadLine();
                        SW.WriteLine(t);
                    }
                    SW.Close();
                }
                SR.Close();
            }
            FI.MoveTo(fullName + ".SKIP");
        }
    }

    public class Ref_Index
    {
        public Dictionary<string, List<Seq>> MainHash;
        public Dictionary<string, Dictionary<short, int>> RefIndexHash;
        public int BitSize;
        public int TotalSequence_Count;
        public int SmallestLibraryCount;
        public int LargestLibraryCount;

        public Ref_Index(Seq_List sList, int BitSize)
        {
            this.BitSize = BitSize;
            MainHash = new Dictionary<string, List<Seq>>();
            RefIndexHash = new Dictionary<string, Dictionary<short, int>>();
            for (int i = 0; i < sList.Pieces.Count; i++)
            {
                BreakAdd(sList.Pieces[i], BitSize, sList.Pieces.GetSourceSeq(i));
            }

            Dictionary<short, int> riCounts = new Dictionary<short, int>();
            foreach (Seq s in sList.List)
            {
                if (!riCounts.ContainsKey(s.Reference_Index))
                    riCounts.Add(s.Reference_Index, 0);
                riCounts[s.Reference_Index]++;
            }
            SortedList<int, short> sorted = new SortedList<int, short>();
            foreach (KeyValuePair<short, int> KVP in riCounts)
            {
                if (!sorted.ContainsKey(KVP.Value))
                    sorted.Add(KVP.Value, KVP.Key);
            }
            SmallestLibraryCount = sorted.Keys[0];
            LargestLibraryCount = sorted.Keys.Last();
            TotalSequence_Count = sList.List.Count;
        }

        public void BreakAdd(string seq, int BitSize, List<Seq> SeqSource)
        {
            for (int i = 0; i < Math.Max((seq.Length - BitSize), 1); i++)
            {
                Add(seq.Substring(i, Math.Min(BitSize, seq.Length)), SeqSource);
            }
        }

        public void Add(string SeqBit, List<Seq> SeqSource)
        {
            if (!MainHash.ContainsKey(SeqBit))
            {
                MainHash.Add(SeqBit, new List<Seq>());
                RefIndexHash.Add(SeqBit, new Dictionary<short, int>());
            }
            if (MainHash[SeqBit].Count < LimitInList)
            {
                MainHash[SeqBit].AddRange(SeqSource);
                foreach (Seq s in SeqSource)
                {
                    if (!RefIndexHash[SeqBit].ContainsKey(s.Reference_Index))
                        RefIndexHash[SeqBit].Add(s.Reference_Index, 0);
                    RefIndexHash[SeqBit][s.Reference_Index]++;
                }
            }
            else
            {

            }
        }

        public const int LimitInList = 50000;

        public void UpdateRefIndexTrack(string Seqbit, Dictionary<short, int> RefIndexTrack)
        {
            if (!RefIndexHash.ContainsKey(Seqbit)) return;
            Dictionary<short, int> tRI_Counts = RefIndexHash[Seqbit];
            foreach (KeyValuePair<short, int> KVP in tRI_Counts)
            {

            }
        }
    }

    public class Seq_List
    {
        public Seq_Schema Schema;
        public Seq_Piece Pieces;
        public List<Seq> List;
        public int MeanAlignSize { get => (int)(SumAlignSize / List.Count); }
        private double SumAlignSize = 0;
        internal static bool SkipRankings = false;

        public HashSet<string> Libraries;

        //public Dictionary<Seq, int> SeqHash;

        public Seq_List()
        {
            Schema = new Seq_Schema();
            Pieces = new Seq_Piece();
            List = new List<Seq>();
            Libraries = new HashSet<string>();
        }

        public void AddSeq(string Name, string seq, short ReferenceIndex, string[] metadata = null, string LibraryName = "")
        {
            if (seq == "GRNA" || seq == "SEQ") return;
            bool SeqRevComp = Library_Aligner.S.ReverseComp_IncomingRefSequence;
            Seq S = new Seq(SeqRevComp ? SeqHelper.ReverseComplement(seq) : seq, Schema, Pieces, ReferenceIndex);
            S.Name = Name;
            S.Library = LibraryName;
            S.ExtraMetaData = string.Join(";", metadata.Skip(2));
            List.Add(S);
            SumAlignSize += S.OriginalSeq.Length;
        }

        public void ImportFromFile(string FilePath, Library_Preference Prefs, StreamWriter SW)
        {
            string[] split = new string[1] { "\t" };
            string t, h, more;
            string[] arr;
            FileInfo FI = new FileInfo(FilePath);
            string LibraryName = Path.GetFileNameWithoutExtension(FI.Name).Trim();
            Libraries.Add(LibraryName);
            short RefIdx = Prefs.Get_orAdd_ReferenceIndex(FI.Name);
            StreamReader SR = new StreamReader(FilePath);
            t = SR.ReadLine();
            bool Continue = true; short d = 0;
            while (Continue)
            {
                if (t.StartsWith(">"))
                {
                    h = t;
                    t = SR.ReadLine(); more = "";
                    while (!t.StartsWith(">"))
                    {
                        more += t;
                        if (SR.EndOfStream) break;
                        if (!SR.EndOfStream) t = SR.ReadLine();
                        else Continue = false;
                    }
                    AddSeq(h, more.ToUpper(), RefIdx);
                    if (SW != null)
                    {
                        SW.WriteLine(FilePath + "\t" + h + "\t" + more.ToUpper());
                    }
                }
                else
                {
                    arr = t.Split(split, StringSplitOptions.None);
                    if (arr.Length >= 2)
                    {
                        if (arr[0].Trim() != "" && !arr[1].Contains("S") && !arr[1].Contains("E")) //This is to ensure that it really is sequence data
                        {
                            if (arr[0].Contains("_"))
                            {
                                if (!short.TryParse(arr[0].Substring(0, arr[0].IndexOf("_")), out d)) d = RefIdx;
                            }
                            AddSeq(arr[0], arr[1].ToUpper(), d, arr, LibraryName);
                        }
                    }
                    if (!SR.EndOfStream) t = SR.ReadLine();
                    else Continue = false;
                }
            }
            SR.Close();
        }

        public string ImportAllFilesInFolder(string FolderPath, Library_Preference LibraryPref = null, string ExportPath_Table = "")
        {
            if (LibraryPref == null) LibraryPref = new Library_Preference();
            try
            {
                Console.WriteLine("Importing Sequences..");
                StreamWriter SW = null;
                if (ExportPath_Table != "") SW = new StreamWriter(ExportPath_Table);
                DirectoryInfo DI = new DirectoryInfo(FolderPath);
                foreach (FileInfo FI in DI.GetFiles("*.txt"))
                {
                    if (!FI.Name.StartsWith("res_") && !FI.Name.StartsWith(".") && !FI.Name.StartsWith("ref_") && !FI.Name.StartsWith("Report"))
                        ImportFromFile(FI.FullName, LibraryPref, SW);
                }
                if (ExportPath_Table != "") SW.Close();
                return "";
            }
            catch (Exception E)
            {
                return "Error importing reference library: " + E.Message;
            }
        }

        public void Export_RefWithResults(Aligner A, string FullName, Library_Preference LibraryPref, string Postfix = "")
        {
            StreamWriter SW = new StreamWriter(FullName);

            List<string> IndicesToReport = LibraryPref.IndexName_StartsWith;
            IndicesToReport.AddRange(LibraryPref.UnMatched);
            SW.WriteLine(Seq.ResultsForExport_Headers(A, IndicesToReport, false, true));
            SW.WriteLine(Seq.ResultsForExport_Headers(A, IndicesToReport, true, false));
            SW.WriteLine(Seq.ResultsForExport_Headers(A, IndicesToReport, false, false));
            for (int i = 0; i < List.Count; i++)
            {
                SW.WriteLine(List[i].ResultsForExport(IndicesToReport) + "\t" + Postfix);
            }
            SW.Close();
        }

        internal void Export_RefLong(Aligner A, string FullName, Library_Preference LibraryPref)
        {
            StreamWriter SW = new StreamWriter(FullName);

            List<string> IndicesToReport = LibraryPref.UnMatched;
            List[0].ResultsExportLong(SW, A, null, IndicesToReport, true);
            for (int i = 0; i < List.Count; i++)
            {
                List[i].ResultsExportLong(SW, A, List[i], IndicesToReport, false);
            }
            SW.Close();
        }

        internal void CalculateRanks(Aligner a, Library_Preference lP)
        {
            Console.WriteLine("Calculating Ranks . .");
            Dictionary<string, SortedList<double, Seq>> RankedList = new Dictionary<string, SortedList<double, Seq>>();
            string SetKey = "";
            double ValKey = 0;
            Random R = new Random();
            foreach (string index in lP.UnMatched)
            {
                for (int i = 0; i < List.Count; i++)
                {
                    SetKey = List[i].Reference_Index + " " + index;
                    if (!RankedList.ContainsKey(SetKey)) RankedList.Add(SetKey, new SortedList<double, Seq>());
                    if (List[i].ResultsSource_Reads.ContainsKey(index))
                    {
                        ValKey = List[i].ResultsSource_Reads[index];
                    }
                    else
                    {
                        ValKey = 0;
                    }
                    while (RankedList[SetKey].ContainsKey(ValKey)) ValKey += R.NextDouble() * 0.000001;
                    RankedList[SetKey].Add(ValKey, List[i]);
                }
                for (int i = 0; i < List.Count; i++)
                {
                    SetKey = List[i].Reference_Index + " " + index;
                    //if (RankedList[SetKey].ContainsValue(List[i])) {
                    List[i].ResultsSource_RankPercent.Add(SetKey, 1 - ((double)RankedList[SetKey].IndexOfValue(List[i]) / RankedList[SetKey].Count));
                }
            }
        }
    }

    public class Seq_Schema
    {
        public List<string> ListStatics;
        public List<int> ListFull;

        public int Count { get => ListFull.Count; }

        public Seq_Schema()
        {
            ListStatics = new List<string>(2);
            ListFull = new List<int>(5);
        }
        public void AddUnique()
        {
            ListFull.Add(eUnique);
        }

        public void AddRedundantBlock()
        {
            ListFull.Add(eRedundantBlock);
        }

        public void AddStatic(string seq)
        {
            ListStatics.Add(seq);
            ListFull.Add(ListStatics.Count - 1);
        }

        public const int eRedundantBlock = -1;
        public const int eUnique = -2;
    }

    public class Seq
    {
        private List<int> Parts;
        private List<string> Uniques;
        public string Name;
        public string Library;
        public bool RepresentsAll;
        public string ExtraMetaData;
        public string OriginalSeq;
        public short Reference_Index;
        public string FullName => Library + "." + Name;
        public string Prefix => Library + "\t" + Reference_Index + "\t" + Name + "\t";
        //public short SubPool;
        public Dictionary<short, double> Result_Reads;
        public Dictionary<string, double> ResultsSource_Reads;
        public Dictionary<string, double> ResultsSource_RankPercent;
        private Object thisLock = new Object();
        //public int Accessed = 0;

        public override string ToString()
        {
            return Name;
        }

        public void AddResults(short ReferenceIndex_Read, double Adjusted_Reads, string Source)
        {
            lock (thisLock)
            {
                if (Result_Reads == null) Result_Reads = new Dictionary<short, double>(2);
                if (Result_Reads.ContainsKey(ReferenceIndex_Read))
                    Result_Reads[ReferenceIndex_Read] += Adjusted_Reads;
                else
                    Result_Reads.Add(ReferenceIndex_Read, Adjusted_Reads);

                if (ResultsSource_Reads == null) ResultsSource_Reads = new Dictionary<string, double>();
                if (ResultsSource_Reads.ContainsKey(Source))
                    ResultsSource_Reads[Source] += Adjusted_Reads;
                else
                    ResultsSource_Reads.Add(Source, Adjusted_Reads);
            }
            //Interlocked.Add(ref Accessed, 1);
        }

        public string ResultsForExport(List<string> IndicesToReport)
        {
            double Reads_Total = 0;
            double Fraction_Correct = 0;
            double Fraction_Second = 0;
            double Reads_Correct = 0;
            string Identity_Second = "";
            double Reads_Second = 0;
            if (Result_Reads != null)
            {
                if (Result_Reads.ContainsKey(Reference_Index))
                {
                    Reads_Correct = Result_Reads[Reference_Index];
                    Reads_Total += Reads_Correct;
                }
                if ((Result_Reads.Count > 1) || (Reads_Total == 0))
                {
                    //First build the sorted list of results with the highest reads first
                    SortedList<double, short> sorted = new SortedList<double, short>(Result_Reads.Count);
                    double tReads;
                    foreach (KeyValuePair<short, double> KVP in Result_Reads)
                    {
                        if (KVP.Key != Reference_Index) //Only put it in here if it is not the main one we are looking for
                        {
                            tReads = KVP.Value;
                            Reads_Total += tReads;
                            while (sorted.ContainsKey(-tReads)) tReads -= 0.000001;
                            sorted.Add(-tReads, KVP.Key);
                        }
                    }
                    //Now put the results for the next best Ref_index
                    if (sorted.Count >= 1)
                    {
                        Identity_Second = sorted.Values[0].ToString();
                        Reads_Second = -sorted.Keys[0];
                        Fraction_Second = Reads_Second / Reads_Total;
                    }
                }
            }
            //if ((Accessed > 0) && (Reads_Total < 0.2)) //This implies a parallelization error
            if (Reads_Total == 0)
                Fraction_Correct = 0;
            else
                Fraction_Correct = Reads_Correct / Reads_Total;
            //Now put together all findings from different indices
            StringBuilder sb = new StringBuilder();
            if (ResultsSource_Reads == null) ResultsSource_Reads = new Dictionary<string, double>();
            for (int i = 0; i < IndicesToReport.Count; i++)
            {
                if (ResultsSource_Reads.ContainsKey(IndicesToReport[i]))
                    sb.Append(ResultsSource_Reads[IndicesToReport[i]] + "\t");
                else
                    sb.Append("0\t");
            }
            string t = Reference_Index + "\t" + Name + "\t" + Fraction_Correct + "\t" + Reads_Correct + "\t" + Identity_Second + "\t" + Fraction_Second + "\t" + Reads_Second + "\t" + sb;
            return t;
        }

        public static string ResultsForExport_Headers(Aligner A, List<string> IndicesToReport, bool Counts, bool Failures)
        {
            StringBuilder sb = new StringBuilder(); string t; string key;
            for (int i = 0; i < IndicesToReport.Count; i++)
            {
                key = IndicesToReport[i];
                if (!A.SourceStats.ContainsKey(key)) t = "";
                else
                {
                    if (Counts) t = A.SourceStats[key].Sequences.ToString();
                    else if (Failures) t = A.SourceStats[key].Sequences_Failed.ToString();
                    else t = key;
                }
                sb.Append(t + "\t");
            }
            string Pre1 = "Reference_Index" + "\t" + "Name" + "\t" + "Fraction_Correct" + "\t" + "Reads_Correct" + "\t" + "Identity_Second" + "\t" + "Fraction_Second" + "\t" + "Reads_Second";
            string Pre2 = "" + "\t" + "" + "\t" + "" + "\t" + "" + "\t" + "" + "\t" + "" + "\t" + (Counts ? "Counts" : "Failures");
            string Pre = (Counts || Failures) ? Pre2 : Pre1;
            return Pre + "\t" + sb;
        }

        internal void ResultsExportLong(StreamWriter SW, Aligner A, Seq RefSeq, List<string> IndicesToReport, bool ExportHeaders)
        {
            //if (RefSeq == null) { RefSeq = new Seq(); RefSeq.Reference_Index = 0; RefSeq.Name = ""; RefSeq.Library = ""; }
            if (ExportHeaders)
            {
                SW.WriteLine("rLibrary" + "\t" + "rSubPool" + "\t" + "rName" + "\t" + "rFullName" + "\t" + "rMeta" + "\t" + PlateWellMeta("") + "\t" + "Reads" + "\t" + "%Reads" + "\t" + "LogNormReads" + "\t" + "%Reads of Mapped" + "\t" + "Mapped Reads" + "\t" + "TotalReads" + "\t" + "%Rank");
                SW.WriteLine("L0" + "\t" + "0" + "\t" + "_Fake" + "\t" + "_Fake" + "\t" + "NA" + "\t" + PlateWellMeta("1-A01-Nada") + "\t" + "0.1" + "\t" + "0.1" + "\t" + "0.1" + "\t" + "0.1" + "\t" + "0.1" + "\t" + "0.1" + "\t" + "0.1");
                if (ResultsSource_Reads == null) ResultsSource_Reads = new Dictionary<string, double>();
                return;
            }
            string Start = RefSeq.Prefix + RefSeq.FullName + "\t" + RefSeq.ExtraMetaData + "\t"; //This is all the Reference-based data, PlateWell Meta contains most of the Index-based-data
            for (int i = 0; i < IndicesToReport.Count; i++)
            {
                double Reads = 0; double LogNormReads = 0;
                double TotalReads = (A.SourceStats[IndicesToReport[i]].Sequences);
                double Failed = (A.SourceStats[IndicesToReport[i]].Sequences_Failed);
                double Mapped = (A.SourceStats[IndicesToReport[i]].Sequences_Mapped);
                double RankPercent = Seq_List.SkipRankings ? -1 : ResultsSource_RankPercent[Reference_Index + " " + IndicesToReport[i]];
                if (ResultsSource_Reads.ContainsKey(IndicesToReport[i]))
                {
                    Reads = ResultsSource_Reads[IndicesToReport[i]];
                    LogNormReads = Math.Log(1 + (10000 * Reads / TotalReads), 2);
                }
                SW.WriteLine(Start + PlateWellMeta(IndicesToReport[i]) + "\t" + Reads + "\t" + (Reads / TotalReads) + "\t" + LogNormReads + "\t" + (Reads / Mapped) + "\t" + Mapped + "\t" + TotalReads + "\t" + RankPercent);
            }
        }

        private static string PlateWellMeta(string IndexName)
        {
            if (IndexName == "")
            {
                return "Index" + "\t" + "iPlate" + "\t" + "iWell" + "\t" + "iWellProper" + "\t" + "iRow" + "\t" + "iCol";
            }
            else
            {
                string[] arr;
                arr = IndexName.Replace("_", "-").Split('-');
                string Row = arr.Length == 1 ? "" : arr[1].Length == 3 ? arr[1].Substring(0, 1) : "";
                string Col = arr.Length == 1 ? "0" : arr[1].Length == 3 ? arr[1].Substring(1, 2) : "0";
                int nCol = 0;
                int.TryParse(Col, out nCol);
                return IndexName + "\t" + arr[0] + "\t" + Row + nCol + "\t" + Row + nCol.ToString("00") + "\t" + Row + "\t" + nCol;
            }
        }

        public Seq()
        {
            RepresentsAll = true;
            Name = "!All";
            OriginalSeq = "";
            Reference_Index = -1;
        }

        public Seq(string seq, Seq_Schema Schema, Seq_Piece sPieces, short ReferenceIndex)
        {
            ResultsSource_RankPercent = new Dictionary<string, double>();
            ResultsSource_Reads = new Dictionary<string, double>();

            seq = seq.Replace(".", ""); seq = seq.Replace("-", "");
            OriginalSeq = seq;
            Reference_Index = ReferenceIndex;
            Parts = new List<int>();
            Uniques = new List<string>();
            RepresentsAll = false;

            int l_start = 0;
            int l_end = 0;
            int l_next_start = 0;
            string t;
            for (int i = 0; i < Schema.Count + 1; i++)
            {
                if (i < Schema.Count)
                {
                    switch (Schema.ListFull[i])
                    {
                        case Seq_Schema.eRedundantBlock:
                        case Seq_Schema.eUnique:
                            l_start = l_next_start;
                            break;
                        default:
                            t = Schema.ListStatics[Schema.ListFull[i]];
                            l_end = seq.IndexOf(t);
                            if (l_end < 0)
                            {
                                l_end = seq.Length; l_next_start = l_end;
                                if (Schema.ListFull[i - 1] == Seq_Schema.eRedundantBlock)
                                {
                                    t = seq.Substring(l_start, l_end - l_start);
                                    Uniques.Add(t);
                                    Parts.Add(sPieces.AddOrGetIndex(t, this));
                                    return;
                                }
                            }
                            else
                            {
                                l_next_start = l_end + t.Length + 1;
                            }
                            break;
                    }
                }
                else
                {
                    l_end = seq.Length;
                }
                if (i > 0)
                {
                    switch (Schema.ListFull[i - 1])
                    {
                        case Seq_Schema.eUnique:
                            t = seq.Substring(l_start, l_end - l_start);
                            Uniques.Add(t);
                            Parts.Add(sPieces.AddOrGetIndex(t, this));
                            break;
                        case Seq_Schema.eRedundantBlock:
                            t = seq.Substring(l_start, l_end - l_start);
                            Parts.Add(sPieces.AddOrGetIndex(t, this));
                            break;
                        default:
                            //Since this is in all, don't add a specific sequence
                            Parts.Add(sPieces.AddOrGetIndex_InAllSeq(Schema.ListStatics[Schema.ListFull[i - 1]]));
                            break;
                    }
                }
                if (l_next_start >= seq.Length) break;
            }
        }
    }

    public class Seq_Piece
    {
        private Dictionary<string, int> PieceHash;
        private List<string> Pieces;
        private List<List<Seq>> SourceList;

        public List<Seq> GetSourceSeq(int index)
        {
            return SourceList[index];
        }

        public int Count { get => Pieces.Count; }
        public string this[int index] { get => Pieces[index]; }

        public Seq_Piece()
        {
            PieceHash = new Dictionary<string, int>();
            Pieces = new List<string>();
            SourceList = new List<List<Seq>>();
        }

        /// <summary>
        /// This one is called if the 
        /// </summary>
        public int AddOrGetIndex_InAllSeq(string SeqPiece)
        {
            int idx;
            if (PieceHash.ContainsKey(SeqPiece))
                idx = PieceHash[SeqPiece];
            else
            {
                Pieces.Add(SeqPiece);
                SourceList.Add(new List<Seq>(1) { new Seq() });
                idx = Pieces.Count - 1;
                PieceHash.Add(SeqPiece, idx);
            }
            return idx;
        }

        public int AddOrGetIndex(string SeqPiece, Seq SeqSource)
        {
            int idx;
            if (PieceHash.ContainsKey(SeqPiece))
                idx = PieceHash[SeqPiece];
            else
            {
                Pieces.Add(SeqPiece);
                SourceList.Add(new List<Seq>());
                idx = Pieces.Count - 1;
                PieceHash.Add(SeqPiece, idx);
            }
            SourceList[idx].Add(SeqSource);
            return idx;
        }
    }
}