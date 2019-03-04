using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using LumenWorks.Framework.IO.Csv;
using MathNet.Numerics.Statistics;
using CSMSL.Proteomics;
using System.IO;


namespace GlycoCompiler
{
    class glycoCompiler
    {
        List<string> inputFilePaths;
        string uniprotGlycoDBFile;
        string outputFolder;
        List<double> glycanMasses = new List<double>();
        double scoreFilter;
        double deltaModFilter;
        double logProbFilter;
        double pepLengthFilter;
        int modCountFilter;

        public glycoCompiler( List<string> inputFilePaths, string uniprotGlycoDBFile, string outputFolder, double scoreFilter,
                                double deltaModFilter, double logProbFilter, double pepLengthFilter, int modCountFilter)
        {
            this.inputFilePaths = inputFilePaths;
            this.uniprotGlycoDBFile = uniprotGlycoDBFile;
            this.outputFolder = outputFolder;
            this.scoreFilter = scoreFilter;
            this.deltaModFilter = deltaModFilter;
            this.logProbFilter = logProbFilter;
            this.pepLengthFilter = pepLengthFilter;
            this.modCountFilter = modCountFilter;
        }

        public void compile()
        {
            OnUpdateProgress(0.1);
            StreamWriter writer = new StreamWriter(@outputFolder + "\\combinedList_Filtered.txt");

            //Concatenate Files
            bool flag = true;
            foreach(string file in inputFilePaths)
            {
                using (CsvReader reader = new CsvReader(new StreamReader(@file), true, '\t'))
                {

                    List<string> headers = reader.GetFieldHeaders().ToList();

                    if (flag)
                    {      
                        string firstLine = "";

                        foreach(string header in headers)
                        {
                            firstLine += header + "\t";
                        }

                        writer.WriteLine(firstLine + "File");

                        flag = false;
                    }

                    while (reader.ReadNextRecord())
                    {
                        string nextLine = "";

                        if (Double.Parse(reader["Delta Mod.Score"]) >= deltaModFilter &&
                                Double.Parse(reader["Score"]) >= scoreFilter &&
                                Double.Parse(reader["|Log Prob|"]) >= logProbFilter &&
                                Double.Parse(reader["Delta Mod.Score"]) >= deltaModFilter &&
                                reader["Peptide"].Length >= pepLengthFilter &&
                                reader["Glycans"].Split(';').Count() <= modCountFilter &&
                                Double.Parse(reader["FDR1D"]) <= 0.01 &&
                                !reader["ProteinName"].Contains("Reverse")
                                )
                        {
                            foreach (string header in headers)
                            {


                                nextLine += reader[header] + "\t";

                            }
                            writer.WriteLine(nextLine + Path.GetFileName(file));
                        }      
                    }
                }      
            }

            writer.Close();

            OnUpdateProgress(0.25);

            // Initialize readers for combined list and uniprot glyco File 
            StreamReader inputFile = new StreamReader(@outputFolder + "\\combinedList_Filtered.txt");
            StreamReader UniprotGlycoDBfile = new StreamReader(@uniprotGlycoDBFile);
            string outputPath = @outputFolder;            

            // Initialize writers
            StreamWriter outputSummary = new StreamWriter(outputPath + "\\_dataSummary_.txt");
            StreamWriter outputPSMs = new StreamWriter(outputPath + "\\_GlycoPSMs.txt");
            StreamWriter outputPeptides = new StreamWriter(outputPath + "\\_GlycoPeptides.txt");
            StreamWriter outputSequences = new StreamWriter(outputPath + "\\_GlycoUniqSeq.txt");
            StreamWriter outputProteins = new StreamWriter(outputPath + "\\_GlycoProteins.txt");
            StreamWriter outputSites = new StreamWriter(outputPath + "\\_GlycoSites.txt");
            StreamWriter outputGlycans = new StreamWriter(outputPath + "\\_Glycans.txt");

            // Read in uniprot glyco data 
            List<UniProtGlycoSite> UniprotGlycoBD = new List<UniProtGlycoSite>();
            using (var csv = new CsvReader(UniprotGlycoDBfile, true))
            {
                while (csv.ReadNextRecord())
                {
                    string uniprotID = csv["UniprotID"];
                    string defLine = csv["DefLine"];
                    double proteinMW = double.Parse(csv["ProteinMW"]);
                    int glycoSite = int.Parse(csv["GlycoSite"]);
                    string glycoType = csv["GlycoType"];
                    string evidenceType = csv["EvidenceType"];
                    int evidenceNumber = int.Parse(csv["EvidenceNumber"]);

                    UniProtGlycoSite databaseEntry = new UniProtGlycoSite(uniprotID, defLine, proteinMW, glycoSite, glycoType, evidenceType, evidenceNumber);
                    UniprotGlycoBD.Add(databaseEntry);
                }
            }

            // Lists for output files
            List<GlycoPSM> glycoPSMs = new List<GlycoPSM>();                                                                    //will include targets and decoys, use for FDR calculation
            List<GlycoPSM> hcdPSMs = new List<GlycoPSM>();                                                                      //will include targets and decoys, use for FDR calculation
            List<GlycoPSM> etdPSMs = new List<GlycoPSM>();                                                                      //will include targets and decoys, use for FDR calculation

            List<GlycoPSM> glycoPSMsNoFilter = new List<GlycoPSM>();                                                            //will include targets and decoys, use for FDR calculation, forgoes 1% FDR Filter at PSM level
            List<GlycoPSM> hcdPSMsNoFilter = new List<GlycoPSM>();                                                              //will include targets and decoys, use for FDR calculation, forgoes 1% FDR Filter at PSM level
            List<GlycoPSM> etdPSMsNoFilter = new List<GlycoPSM>();                                                              //will include targets and decoys, use for FDR calculation, forgoes 1% FDR Filter at PSM level

            List<GlycoPSM> glycoPSMsLocalized = new List<GlycoPSM>();
            List<GlycoPSM> hcdPSMsLocalized = new List<GlycoPSM>();
            List<GlycoPSM> etdPSMsLocalized = new List<GlycoPSM>();

            Dictionary<string, List<GlycoPSM>> uniqueSequencesDictionary = new Dictionary<string, List<GlycoPSM>>();            //string is peptide sequence only

            // Peptides to be parsed for highest scoring PSM
            Dictionary<string, List<GlycoPSM>> glycoPeptideDictionary = new Dictionary<string, List<GlycoPSM>>();               //string is peptides sequence with mass mods at locations
            Dictionary<string, List<GlycoPSM>> glycoSitesDictionary = new Dictionary<string, List<GlycoPSM>>();                 //string is protein UniProtID concatenated with site
            Dictionary<string, List<GlycoPSM>> glycansDictionary = new Dictionary<string, List<GlycoPSM>>();                    //string is glycan, i.e., HexNAc(4)Hex(3), etc.
            Dictionary<string, List<GlycoPSM>> glycoProteinsDictionary = new Dictionary<string, List<GlycoPSM>>();              //string is UniprotID
            Dictionary<string, HashSet<int>> glycoProteinsWithSites = new Dictionary<string, HashSet<int>>();                   //string is UniprotID

            using (var csv = new CsvReader(inputFile, true, '\t'))
            {
                while (csv.ReadNextRecord())
                {
                    List<string> headers = csv.GetFieldHeaders().ToList();

                    List<int> glycanPositionsList = new List<int>();
                    List<string> mods = new List<string>();
                    List<string> glycans = new List<string>();

                    /*
                    outputList.Write("PID\tProt.Rank\tSequence\tPeptideParseFriendly\tPos.\tMods(variable)\tGlycans\tPEP2D\tPEP1D\t|Log Prob|\tScore\tDeltaScore\tDelta Mod.Score\tz\tObs.m/z\t");
                    outputList.Write("Calc.m/z\tppmerr.\tObs.m/z\tCalc.m/z\tCleavage\tGlycansPos.\tProteinName\tProt.Id\tScanTime\tScan #\tMods(fixed)\tFDR2D\tFDR1D\t");
                    outputList.WriteLine("FDR uniq.2D\tFDR uniq.1D\tq-value2D\tq-value1D\tisGlycoPeptide\tmodsPassedCheck\tpositionPassedCheck\tfragmentation\tMasterScan\tPeak126\tPeak138\tPeak144\tPeak168\tPeak186\tPeak204\tPeak274\tPeak292\tPeak366\tGlcNAc/GalNAcratio");
                    */

                    string PID = "NA";   //csv["PID"];
                    string protRank = csv["Prot.Rank"];
                    string sequence = csv["Sequence"];
                    string peptidesToBeParsed = csv["PeptideParseFriendly"];
                    int peptideStartPosition = int.Parse(csv["Pos."]);
                    string modsToBeParsed = csv["Mods(variable)"];
                    string glycansToBeParsed = csv["Glycans"];
                    double PEP2D = double.Parse(csv["PEP2D"]);
                    double PEP1D = double.Parse(csv["PEP1D"]);
                    double logProb = Math.Abs(Math.Log10(PEP1D));
                    double score = double.Parse(csv["Score"]);
                    double deltaScore = double.Parse(csv["DeltaScore"]);
                    double deltaModScore = double.Parse(csv["Delta Mod.Score"]);
                    int charge = int.Parse(csv["Charge"]);
                    double mzObs = double.Parse(csv["Obs.m/z"]);
                    double mzCalc = double.Parse(csv["Precursor Theoretical m/z (Th)"]);
                    double ppmError = double.Parse(csv["ppmerr."]);
                    double obsMH = double.Parse(csv["Obs.MH"]);
                    double calcMH = double.Parse(csv["Calc.MH"]);
                    string cleavage = csv["Cleavage"];
                    string glycanPositions = csv["GlycansPos."];
                    string proteinName = csv["ProteinName"];
                    int protID = int.Parse(csv["Prot.Id"]);
                    double scanTime = double.Parse(csv["ScanTime"]);
                    int scanNumber = int.Parse(csv["Spectrum Number"]);
                    string modsFixed = csv["Mods(fixed)"];
                    double FDR2D = double.Parse(csv["FDR2D"]);
                    double FDR1D = double.Parse(csv["FDR1D"]);
                    double FDR2Dunique = double.Parse(csv["FDR uniq.2D"]);
                    double FDR1Dunique = double.Parse(csv["FDR uniq.1D"]);
                    double qvalue2D = double.Parse(csv["q-value2D"]);
                    double qvalue1D = double.Parse(csv["q-value1D"]);
                    bool isGlycoPeptide = bool.Parse(csv["isGlycoPeptide"]);
                    string fragmentation = csv["fragmentation"];
                    int masterScan = int.Parse(csv["MasterScan"]);
                    string byonicIntensity = "NA"; //csv["ByonicIntensity"];
                    string file = csv["File"];

                    double peak126 = 0;
                    double peak138 = 0;
                    double peak144 = 0;
                    double peak168 = 0;
                    double peak186 = 0;
                    double peak204 = 0;
                    double peak274 = 0;
                    double peak292 = 0;
                    double peak366 = 0;
                    double GlcNAcGalNAcratio = 0;

                    /**
                    double peak126 = double.Parse(csv["Peak126"]);
                    double peak138 = double.Parse(csv["Peak138"]);
                    double peak144 = double.Parse(csv["Peak144"]);
                    double peak168 = double.Parse(csv["Peak168"]);
                    double peak186 = double.Parse(csv["Peak186"]);
                    double peak204 = double.Parse(csv["Peak204"]);
                    double peak274 = double.Parse(csv["Peak274"]);
                    double peak292 = double.Parse(csv["Peak292"]);
                    double peak366 = double.Parse(csv["Peak366"]);                    
                    double GlcNAcGalNAcratio = double.Parse(csv["GlcNAc/GalNAcratio"]);
                    **/

                    double intensity = double.Parse(csv["LFQ Intensity"]);
                    

                    bool seenWithHCD = false;
                    bool seenWithETD = false;
                    bool NXSmotif = false;
                    bool NXTmotif = false;
                    bool isLocalized = false;
                    bool Nlinked = false;
                    bool Olinked = false;

                    bool matchedToUniprot = false;
                    string uniprotEvidenceType = "None";
                    int uniprotEvidenceNumber = 0;

                    string[] parsedPeptide = peptidesToBeParsed.Split(',');
                    Peptide peptide = new Peptide(parsedPeptide[0]);
                    double peptideMonoMass = peptide.MonoisotopicMass;

                    char[] peptideTermini = parsedPeptide[1].ToCharArray();
                    char peptideCterminusNextResidue = peptideTermini[1];

                    string[] parsedSequenceWithMods = sequence.Split('.');
                    string sequenceWithMods = parsedSequenceWithMods[1];

                    string[] parsedProteinName = proteinName.Split('|');
                    string uniprotID = parsedProteinName[1];

                    if (fragmentation.Equals("HCD"))
                        seenWithHCD = true;

                    if (fragmentation.Equals("ETD"))
                        seenWithETD = true;

                    if (!String.IsNullOrEmpty(glycanPositions))
                    {
                        string[] glycansPosParsedArray = glycanPositions.Split(';');
                        for (int i = 0; i < glycansPosParsedArray.Length; i++)
                        {
                            int glycanPos = Convert.ToInt32(glycansPosParsedArray[i]);
                            glycanPositionsList.Add(glycanPos);
                        }
                    }

                    if (!String.IsNullOrEmpty(glycansToBeParsed))
                    {
                        string[] glycansParsedArrary = glycansToBeParsed.Split(';');

                        for (int i = 0; i < glycansParsedArrary.Length; i++)
                        {
                            string glycan = glycansParsedArrary[i];
                            if (glycan[0].Equals(' '))
                            {
                                glycan = glycan.Substring(1);
                            }
                            glycans.Add(glycan);
                        }
                    }

                    int numberOfSites = glycans.Count;



                    if (isGlycoPeptide && !proteinName.Contains("DECOY"))          
                    {

                        if (!String.IsNullOrEmpty(modsToBeParsed))
                        {
                            string[] modsParsedArrary = modsToBeParsed.Split(';');
                            for (int i = 0; i < modsParsedArrary.Length; i++)
                            {
                                string mod = modsParsedArrary[i];
                                if (mod[0].Equals(' '))
                                {
                                    mod = mod.Substring(1);
                                }

                                mods.Add(mod);
                                string modName = GetModName(mod);
                                double modMass = GetModMass(mod);
                                int modPosition = GetModPosition(mod);
                                if (modName.Contains("Glycan"))
                                {
                                    modName = modName + "_" + modMass;
                                }

                                if (modName.Contains("NGlycan"))
                                    Nlinked = true;

                                if (modName.Contains("OGlycan"))
                                    Olinked = true;

                                Modification modToAdd = new Modification(modMass, modName);
                                peptide.AddModification(modToAdd, modPosition);
                                
                                

                                /*
                                Console.WriteLine(mod);
                                Console.WriteLine(modName);
                                Console.WriteLine(modMass);
                                Console.WriteLine(modPosition);
                                Console.WriteLine(peptide.SequenceWithModifications);
                                Console.WriteLine(uniprotEvidenceNumber);
                                */

                                if (modName.Contains("NGlycan") || modName.Contains("OGlycan"))
                                {
                                    glycanMasses.Add(modMass);

                                    if ((modPosition + 2) > peptide.Length)
                                    {
                                        if (!peptide.GetResidue(peptide.Length - 1).Equals('P'))
                                        {
                                            if (peptideCterminusNextResidue.Equals('S'))
                                            {
                                                NXSmotif = true;
                                            }
                                            if (peptideCterminusNextResidue.Equals('T'))
                                            {
                                                NXTmotif = true;
                                            }
                                        }
                                    }
                                    else
                                    {
                                        if (!peptide.GetResidue(modPosition).Letter.Equals('P'))
                                        {
                                            if (peptide.GetResidue(modPosition + 1).Letter.Equals('S'))
                                            {
                                                NXSmotif = true;
                                            }
                                            if (peptide.GetResidue(modPosition + 1).Letter.Equals('T'))
                                            {
                                                NXTmotif = true;
                                            }
                                        }
                                    }
                                }

                                if (modName.Contains("NGlycan") || modName.Contains("OGlycan"))
                                {
                                    if (!glycoProteinsWithSites.ContainsKey(uniprotID))
                                    {
                                        glycoProteinsWithSites.Add(uniprotID, new HashSet<int>());
                                    }
                                    glycoProteinsWithSites[uniprotID].Add(modPosition + peptideStartPosition - 1);

                                    foreach (UniProtGlycoSite site in UniprotGlycoBD)
                                    {
                                        if (site.uniprotID.Equals(uniprotID))
                                        {
                                            if (site.glycoSite == modPosition + peptideStartPosition - 1)
                                            {
                                                matchedToUniprot = true;
                                                uniprotEvidenceType = site.evidenceType;
                                                uniprotEvidenceNumber = site.evidenceNumber;
                                            }
                                        }
                                    }
                                }
                            }
                        }


                        if (!String.IsNullOrEmpty(modsFixed))
                        {
                            string[] modsParsedArrary = modsFixed.Split(';');
                            for (int i = 0; i < modsParsedArrary.Length; i++)
                            {
                                string mod = modsParsedArrary[i];
                                if (mod[0].Equals(' '))
                                {
                                    mod = mod.Substring(1);
                                }
                                mods.Add(mod);
                                string modName = GetModName(mod);
                                double modMass = GetModMass(mod);
                                int modPosition = GetModPosition(mod);

                                Modification modToAdd = new Modification(modMass, modName);
                                peptide.AddModification(modToAdd, modPosition);
                            }
                        }

                        if (deltaModScore >= 10)
                            isLocalized = true;
                        
                        

                        GlycoPSM psm = new GlycoPSM(peptide, peptideMonoMass, peptide.SequenceWithModifications, mods, 
                            glycans, glycanPositionsList, uniprotID, PEP2D, PEP1D, logProb, score, deltaScore, deltaModScore, 
                            mzObs, mzCalc, charge, numberOfSites, ppmError, obsMH, calcMH, cleavage, proteinName, peptideStartPosition, 
                            scanTime, scanNumber, FDR2D, FDR1D, FDR2Dunique, FDR1Dunique, qvalue2D, qvalue1D, fragmentation, masterScan, 
                            peak126, peak138, peak144, peak168, peak186, peak204, peak274, peak292, peak366, GlcNAcGalNAcratio, isGlycoPeptide, 
                            seenWithHCD, seenWithETD, NXSmotif, NXTmotif, isLocalized, Nlinked, Olinked, matchedToUniprot, uniprotEvidenceType, 
                            uniprotEvidenceNumber, intensity, byonicIntensity, file);


                        if (!uniqueSequencesDictionary.ContainsKey(peptide.Sequence))
                        {
                            uniqueSequencesDictionary.Add(peptide.Sequence, new List<GlycoPSM>());
                        }
                        uniqueSequencesDictionary[peptide.Sequence].Add(psm);

                        if (!glycoPeptideDictionary.ContainsKey(peptide.SequenceWithModifications + " " + uniprotID))
                        {
                            glycoPeptideDictionary.Add(peptide.SequenceWithModifications + " " + uniprotID, new List<GlycoPSM>());
                        }
                        glycoPeptideDictionary[peptide.SequenceWithModifications + " " + uniprotID].Add(psm);

                        foreach (int site in glycanPositionsList)
                        {
                            int glycoSite = site + peptideStartPosition - 1;
                            string siteID = uniprotID + "_" + glycoSite;
                            if (!glycoSitesDictionary.ContainsKey(siteID))
                            {
                                glycoSitesDictionary.Add(siteID, new List<GlycoPSM>());
                            }
                            glycoSitesDictionary[siteID].Add(psm);
                            //Console.WriteLine(siteID);
                            //Console.ReadKey();
                        }

                        foreach (string glycan in glycans)
                        {
                            if (!glycansDictionary.ContainsKey(glycan))
                            {
                                glycansDictionary.Add(glycan, new List<GlycoPSM>());
                            }
                            glycansDictionary[glycan].Add(psm);
                            //Console.WriteLine(glycan);
                            //Console.ReadKey();
                        }

                        if (!glycoProteinsDictionary.ContainsKey(uniprotID))
                        {
                            glycoProteinsDictionary.Add(uniprotID, new List<GlycoPSM>());
                            glycoProteinsDictionary[uniprotID].Add(psm);
                        }
                        else
                        {
                            glycoProteinsDictionary[uniprotID].Add(psm);
                        }
                        
                    }

                    if (peptide.Length > 4 && FDR2D <= 0.01 && isGlycoPeptide && score >= scoreFilter)
                    {
                        GlycoPSM psmTargetORDecoy = new GlycoPSM(peptide, peptideMonoMass, peptide.SequenceWithModifications, mods, glycans, 
                                glycanPositionsList, uniprotID, PEP2D, PEP1D, logProb, score, deltaScore, deltaModScore, mzObs, mzCalc, charge, 
                                numberOfSites, ppmError, obsMH, calcMH, cleavage, proteinName, peptideStartPosition, scanTime, scanNumber, FDR2D, 
                                FDR1D, FDR2Dunique, FDR1Dunique, qvalue2D, qvalue1D, fragmentation, masterScan, peak126, peak138, peak144,
                                peak168, peak186, peak204, peak274, peak292, peak366, GlcNAcGalNAcratio, isGlycoPeptide, seenWithHCD, seenWithETD, 
                                NXSmotif, NXTmotif, isLocalized, Nlinked, Olinked, matchedToUniprot, uniprotEvidenceType, uniprotEvidenceNumber, 
                                intensity, byonicIntensity, file);

                        glycoPSMs.Add(psmTargetORDecoy);
                        if (fragmentation.Equals("HCD"))
                            hcdPSMs.Add(psmTargetORDecoy);
                        if (fragmentation.Equals("ETD"))
                            etdPSMs.Add(psmTargetORDecoy);

                        if (isLocalized)
                        {
                            glycoPSMsLocalized.Add(psmTargetORDecoy);
                            if (fragmentation.Equals("HCD"))
                                hcdPSMsLocalized.Add(psmTargetORDecoy);
                            if (fragmentation.Equals("ETD"))
                                etdPSMsLocalized.Add(psmTargetORDecoy);
                        }
                    }

                    GlycoPSM psmTargetORDecoyNoFilter = new GlycoPSM(peptide, peptideMonoMass, peptide.SequenceWithModifications, mods, glycans, 
                        glycanPositionsList, uniprotID, PEP2D, PEP1D, logProb, score, deltaScore, deltaModScore, mzObs, mzCalc, charge, numberOfSites, 
                        ppmError, obsMH, calcMH, cleavage, proteinName, peptideStartPosition, scanTime, scanNumber, FDR2D, FDR1D, FDR2Dunique, FDR1Dunique, 
                        qvalue2D, qvalue1D, fragmentation, masterScan, peak126, peak138, peak144, peak168, peak186, peak204, peak274, peak292, peak366, 
                        GlcNAcGalNAcratio, isGlycoPeptide, seenWithHCD, seenWithETD, NXSmotif, NXTmotif, isLocalized, Nlinked, Olinked, matchedToUniprot, 
                        uniprotEvidenceType, uniprotEvidenceNumber, intensity, byonicIntensity, file);

                    glycoPSMsNoFilter.Add(psmTargetORDecoyNoFilter);
                    if (fragmentation.Equals("HCD"))
                        hcdPSMsNoFilter.Add(psmTargetORDecoyNoFilter);
                    if (fragmentation.Equals("ETD"))
                        etdPSMsNoFilter.Add(psmTargetORDecoyNoFilter);
                }
            }

            OnUpdateProgress(0.5);

            /////////////////////////////////////////////////////////////calculate FDR for all, HCD, and ETD
            double fdrAllGlycoPSMs = CalculateFDR(glycoPSMs);
            double fdrHCD = CalculateFDR(hcdPSMs);
            double fdrETD = CalculateFDR(etdPSMs);

            double fdrAllGlycoPSMsNoFilter = CalculateFDR(glycoPSMsNoFilter);
            double fdrHCDNoFilter = CalculateFDR(hcdPSMsNoFilter);
            double fdrETDNoFilter = CalculateFDR(etdPSMsNoFilter);

            double fdrAllLocalized = CalculateFDR(glycoPSMsLocalized);
            double fdrHCDLocalized = CalculateFDR(hcdPSMsLocalized);
            double fdrETDLocalized = CalculateFDR(etdPSMsLocalized);

            ////////////////////////////////////////////////////////////count and print glycans

            outputGlycans.WriteLine("Glycan\tMass\tLocalized\tUniqueSeq\tUniqueSeqLocalized\tProteins\tProteinsLocalized\tLocalizedUniprotIDs");

            int countOfGlycansAll = glycansDictionary.Count;
            int countOfGlycansLocalized = 0;
            int countOfGlcyansAllNlinked = 0;
            int countOfGlycansLocazliedNlinked = 0;
            int countOfGlycansAllOlinked = 0;
            int countOfGlycansLocalizedOlinked = 0;

            foreach (KeyValuePair<string, List<GlycoPSM>> entry in glycansDictionary)
            {
                int numberOfPSMs = 0;
                int numberOfPSMsLocalized = 0;
                string proteinsLocalized = "NotLocalized";

                HashSet<string> uniqueSequences = new HashSet<string>();
                HashSet<string> uniqueSequencesLocalized = new HashSet<string>();
                HashSet<string> uniprotIDs = new HashSet<string>();
                HashSet<string> uniprotIDsLocalized = new HashSet<string>();

                bool glycanLocalized = false;
                bool Nglycan = false;
                bool NglycanLocalized = false;
                bool Oglycan = false;
                bool OglycanLocalized = false;

                List<double> glycanMasses = new List<double>();

                foreach (GlycoPSM psm in entry.Value)
                {
                    uniqueSequences.Add(psm.peptide.Sequence);
                    uniprotIDs.Add(psm.uniprotID);
                    if (psm.Nlinked)
                        Nglycan = true;
                    if (psm.Olinked)
                        Oglycan = true;

                    if (psm.isLocalized)
                    {
                        glycanLocalized = true;
                        if (psm.Nlinked)
                            NglycanLocalized = true;
                        if (psm.Olinked)
                            OglycanLocalized = true;
                        uniqueSequencesLocalized.Add(psm.peptide.Sequence);
                        uniprotIDsLocalized.Add(psm.uniprotID);

                        /*
                        if (proteinsLocalized.Equals("NotLocalized"))
                        {
                            proteinsLocalized = psm.uniprotID;
                        }
                        else
                        {
                            proteinsLocalized += ";" + psm.uniprotID;
                        }
                        */

                    }


                    foreach (string mod in psm.mods)
                    {
                        string name = GetModName(mod);
                        if (name.Contains("Glycan"))
                        {
                            double modMass = GetModMass(mod);
                            glycanMasses.Add(modMass);
                        }
                    }
                }

                foreach (string protein in uniprotIDsLocalized)
                {
                    if (proteinsLocalized.Equals("NotLocalized"))
                    {
                        proteinsLocalized = protein;
                    }
                    else
                    {
                        proteinsLocalized += ";" + protein;
                    }
                }

                OnUpdateProgress(0.75);

                int numberOfUniqueSequences = uniqueSequences.Count;
                int numberOfUniqueSequencesLocalized = uniqueSequencesLocalized.Count;
                int numberOfProteins = uniprotIDs.Count;
                int numberOfProteinsLocalized = uniprotIDsLocalized.Count;

                var glycanMassesGroup = glycanMasses.GroupBy(v => v);
                int maxCount = glycanMassesGroup.Max(g => g.Count());
                double mode = glycanMassesGroup.First(g => g.Count() == maxCount).Key; // this is the glycan mass

                //Console.WriteLine(proteinsLocalized);
                //Console.ReadKey();


                outputGlycans.WriteLine(entry.Key + "\t" + mode + "\t" + glycanLocalized + "\t" + numberOfUniqueSequences + "\t" + numberOfUniqueSequencesLocalized + "\t" + numberOfProteins + "\t" + numberOfProteinsLocalized + "\t" + proteinsLocalized);


                if (glycanLocalized)
                    countOfGlycansLocalized++;
                if (Nglycan)
                    countOfGlcyansAllNlinked++;
                if (NglycanLocalized)
                    countOfGlycansLocazliedNlinked++;
                if (Oglycan)
                    countOfGlycansAllOlinked++;
                if (OglycanLocalized)
                    countOfGlycansLocalizedOlinked++;
            }

            ////////////////////////////////////////////////////////////count and print glycosites

            outputSites.WriteLine("Protein\tDescription\tSite\tLocalized\tGlycoType\tPSMs\tUniqueSeq\tNo.OfGlycans\tSeenInUniprot\tEvidenceType\tEvidence#\tGlycans");

            int countOfGlycoSitesAll = glycoSitesDictionary.Count;
            int countOfGlycoSitesLocalized = 0;
            int countOfGlycoSitesAllNlinked = 0;
            int countOfGlycoSitesLocazliedNlinked = 0;
            int countOfGlycoSitesAllOlinked = 0;
            int countOfGlycoSitesLocalizedOlinked = 0;

            int countOfSitesSeenInUniprot = 0;
            int countOfSitesSeenInUniprotLocalized = 0;
            int countBySimilarity = 0;
            int countBySimilarityLocalized = 0;
            int countDiscussedNoEvidence = 0;
            int countDiscussedNoEvidenceLocalized = 0;
            int countInferredByCurator = 0;
            int countInferredByCuratorLocalized = 0;
            int countInferredEvidence = 0;
            int countInferredEvidenceLocalized = 0;
            int countReportInLit = 0;
            int countReportInLitLocalized = 0;
            int countSeqAnalysis = 0;
            int countSeqAnalysisLocalized = 0;


            foreach (KeyValuePair<string, List<GlycoPSM>> entry in glycoSitesDictionary)
            {
                string[] parsedGlycoSite = entry.Key.Split('_');

                string uniprotID = parsedGlycoSite[0];
                int site = Convert.ToInt32(parsedGlycoSite[1]);

                bool siteIsLocalized = false;
                bool NlinkedSite = false;
                bool siteIsLocalizedNlinked = false;
                bool OlinkedSite = false;
                bool siteIsLocalizedOlinked = false;

                int numberOfPSMs = 0;
                string glycans = "NoGlycans";
                string description = "needsDescription";
                string glycoType = "blank";

                bool inUniprot = false;
                string evidence = "NotSeenInUniprot";
                int evidenceNumber = 0;

                bool inUniprotLocalized = false;
                bool bySimilarity = false;
                bool bySimilarityLocalized = false;
                bool discussedNoEvidence = false;
                bool discussedNoEvidenceLocalized = false;
                bool inferredByCurator = false;
                bool inferredByCuratorLocalized = false;
                bool inferredByEvidence = false;
                bool inferredByEvidenceLocalized = false;
                bool reportInLit = false;
                bool reportInLitLocalized = false;
                bool seqAnalysis = false;
                bool seqAnalysisLocalized = false;


                HashSet<string> uniqueSequences = new HashSet<string>();
                HashSet<string> glycansAtThisSite = new HashSet<string>();

                foreach (GlycoPSM psm in entry.Value)
                {
                    if (psm.isLocalized)
                    {
                        siteIsLocalized = true;
                        if (psm.Nlinked)
                        {
                            siteIsLocalizedNlinked = true;
                        }
                        else if (psm.Olinked)
                        {
                            siteIsLocalizedOlinked = true;
                        }

                        if (psm.evidenceNumber == 250)
                            bySimilarityLocalized = true;

                        if (psm.evidenceNumber == 303)
                            discussedNoEvidenceLocalized = true;

                        if (psm.evidenceNumber == 305)
                            inferredByCuratorLocalized = true;

                        if (psm.evidenceNumber == 244)
                            inferredByEvidenceLocalized = true;

                        if (psm.evidenceNumber == 269)
                            reportInLitLocalized = true;

                        if (psm.evidenceNumber == 255)
                            seqAnalysisLocalized = true;
                    }
                    numberOfPSMs++;
                    uniqueSequences.Add(psm.peptide.Sequence);

                    foreach (string glycan in psm.glycans)
                    {
                        if (glycans.Equals("NoGlycans"))
                        {
                            glycans = glycan;
                            glycansAtThisSite.Add(glycan);
                        }
                        else
                        {
                            if (!glycansAtThisSite.Contains(glycan))
                            {
                                glycans += ";" + glycan;
                                glycansAtThisSite.Add(glycan);
                            }
                        }
                    }

                    if (psm.matchedToUniprot)
                        inUniprot = true;
                    evidence = psm.evidenceType;
                    evidenceNumber = psm.evidenceNumber;

                    if (psm.evidenceNumber == 250)
                        bySimilarity = true;

                    if (psm.evidenceNumber == 303)
                        discussedNoEvidence = true;

                    if (psm.evidenceNumber == 305)
                        inferredByCurator = true;

                    if (psm.evidenceNumber == 244)
                        inferredByEvidence = true;

                    if (psm.evidenceNumber == 269)
                        reportInLit = true;

                    if (psm.evidenceNumber == 255)
                        seqAnalysis = true;

                    if (psm.Nlinked)
                    {
                        glycoType = "Nlinked";
                        NlinkedSite = true;
                    }
                    else if (psm.Olinked)
                    {
                        glycoType = "Olinked";
                        OlinkedSite = true;
                    }

                    description = psm.proteinName;
                }

                if (siteIsLocalized)
                    countOfGlycoSitesLocalized++;
                if (NlinkedSite)
                    countOfGlycoSitesAllNlinked++;
                if (siteIsLocalizedNlinked)
                    countOfGlycoSitesLocazliedNlinked++;
                if (OlinkedSite)
                    countOfGlycoSitesAllOlinked++;
                if (siteIsLocalizedOlinked)
                    countOfGlycoSitesLocalizedOlinked++;

                if (bySimilarity)
                {
                    countBySimilarity++;
                    countOfSitesSeenInUniprot++;
                }
                if (bySimilarityLocalized)
                {
                    countBySimilarityLocalized++;
                    countOfSitesSeenInUniprotLocalized++;
                }

                if (discussedNoEvidence)
                {
                    countDiscussedNoEvidence++;
                    countOfSitesSeenInUniprot++;
                }
                if (discussedNoEvidenceLocalized)
                {
                    countDiscussedNoEvidenceLocalized++;
                    countOfSitesSeenInUniprotLocalized++;
                }

                if (inferredByCurator)
                {
                    countInferredByCurator++;
                    countOfSitesSeenInUniprot++;
                }
                if (inferredByCuratorLocalized)
                {
                    countInferredByCuratorLocalized++;
                    countOfSitesSeenInUniprotLocalized++;
                }

                if (inferredByEvidence)
                {
                    countInferredEvidence++;
                    countOfSitesSeenInUniprot++;
                }
                if (inferredByEvidenceLocalized)
                {
                    countInferredEvidenceLocalized++;
                    countOfSitesSeenInUniprotLocalized++;
                }

                if (reportInLit)
                {
                    countReportInLit++;
                    countOfSitesSeenInUniprot++;
                }
                if (reportInLitLocalized)
                {
                    countReportInLitLocalized++;
                    countOfSitesSeenInUniprotLocalized++;
                }

                if (seqAnalysis)
                {
                    countSeqAnalysis++;
                    countOfSitesSeenInUniprot++;
                }
                if (seqAnalysisLocalized)
                {
                    countSeqAnalysisLocalized++;
                    countOfSitesSeenInUniprotLocalized++;
                }

                glycans.Count();

                outputSites.WriteLine(uniprotID + "\t" + description + "\t" + site + "\t" + siteIsLocalized + "\t" + glycoType + "\t" + numberOfPSMs + "\t" + uniqueSequences.Count + "\t" + glycansAtThisSite.Count + "\t" + inUniprot + "\t" + evidence + "\t" + evidenceNumber + "\t" + glycans);
            }

            ////////////////////////////////////////////////////////////count and print glycoproteins

            outputProteins.WriteLine("UniprotID\tDescription\tInUniprotAsGlycoProtein\tLocalizedPSMs\tLocalizedUniqPeps\tLocalizedUniqSeq\tLocalizedGlycoSites\tNo.OfLocalizedGlycans\tLocalizedGlycans\t" +
                "AllGlycoPSMs\tAllGlycoUniqPeps\tAllGlycoUniqSeq\tAllGlycoSites\tNo.OfAllGlycans\tAllGlycans");

            int countOfGlycoProteins = glycoProteinsDictionary.Count;
            int countOfGlycoProteinsLocalized = 0;
            int countOfGlycoProteinsNlinked = 0;
            int countOfGlycoProteinsNlinkedLocalized = 0;
            int countOfGlycoProteinsOlinked = 0;
            int countOfGlycoProteinsOlinkedLocalized = 0;

            foreach (KeyValuePair<string, List<GlycoPSM>> entry in glycoProteinsDictionary)
            {
                string uniprotID = entry.Key;
                string description = "blank";

                int numberOfPSMs = 0;
                int numberOfPSMsLocalized = 0;

                HashSet<string> uniquePeptides = new HashSet<string>();
                HashSet<string> uniqueSequences = new HashSet<string>();
                HashSet<string> glycans = new HashSet<string>();
                HashSet<int> glycoSites = new HashSet<int>();

                HashSet<string> uniquePeptidesLocalized = new HashSet<string>();
                HashSet<string> uniqueSequencesLocalized = new HashSet<string>();
                HashSet<string> glycansLocalized = new HashSet<string>();
                HashSet<int> glycoSitesLocalized = new HashSet<int>();

                bool inUniprotAsGlycoProtein = false;

                bool localized = false;
                bool Nlinked = false;
                bool NlinkedLocalized = false;
                bool Olinked = false;
                bool OlinkedLocalized = false;

                string glycansSeenOnThisProtein = "none";
                string glycansSeenOnThisProteinLocalized = "none";

                foreach (GlycoPSM psm in entry.Value)
                {
                    if (psm.matchedToUniprot)
                        inUniprotAsGlycoProtein = true;

                    description = psm.proteinName;

                    if (psm.isLocalized)
                    {
                        numberOfPSMsLocalized++;
                        localized = true;
                        uniquePeptidesLocalized.Add(psm.sequenceWithMods);
                        uniqueSequencesLocalized.Add(psm.peptide.Sequence);
                        foreach (int site in psm.glycanPositions)
                        {
                            int actualPosition = site + psm.peptideStartPosition - 1;
                            glycoSitesLocalized.Add(actualPosition);
                        }
                        foreach (string glycan in psm.glycans)
                        {
                            glycansLocalized.Add(glycan);
                        }

                        if (psm.Nlinked)
                            NlinkedLocalized = true;
                        if (psm.Olinked)
                            OlinkedLocalized = true;
                    }

                    if (psm.Nlinked)
                        Nlinked = true;
                    if (psm.Olinked)
                        Olinked = true;

                    numberOfPSMs++;
                    uniquePeptides.Add(psm.sequenceWithMods);
                    uniqueSequences.Add(psm.peptide.Sequence);
                    foreach (int site in psm.glycanPositions)
                    {
                        int actualPosition = site + psm.peptideStartPosition - 1;
                        glycoSites.Add(actualPosition);
                    }
                    foreach (string glycan in psm.glycans)
                    {
                        glycans.Add(glycan);
                    }
                }

                if (localized)
                    countOfGlycoProteinsLocalized++;
                if (Nlinked)
                    countOfGlycoProteinsNlinked++;
                if (NlinkedLocalized)
                    countOfGlycoProteinsNlinkedLocalized++;
                if (Olinked)
                    countOfGlycoProteinsOlinked++;
                if (OlinkedLocalized)
                    countOfGlycoProteinsOlinkedLocalized++;

                foreach (string glycan in glycans)
                {
                    if (glycansSeenOnThisProtein.Equals("none"))
                        glycansSeenOnThisProtein = glycan;
                    else
                        glycansSeenOnThisProtein += ";" + glycan;
                }

                foreach (string glycanLocalized in glycansLocalized)
                {
                    if (glycansSeenOnThisProteinLocalized.Equals("none"))
                        glycansSeenOnThisProteinLocalized = glycanLocalized;
                    else
                        glycansSeenOnThisProteinLocalized += ";" + glycanLocalized;
                }

                if (description.Equals("blank"))
                {
                    int asdf = 0;
                }

                outputProteins.WriteLine(uniprotID + "\t" + description + "\t" + inUniprotAsGlycoProtein + "\t" + numberOfPSMsLocalized + "\t" + uniquePeptidesLocalized.Count + "\t" + uniqueSequencesLocalized.Count + "\t" + glycoSitesLocalized.Count + "\t" +
                    glycansLocalized.Count + "\t" + glycansSeenOnThisProteinLocalized + "\t" + numberOfPSMs + "\t" + uniquePeptides.Count + "\t" + uniqueSequences.Count + "\t" + glycoSites.Count + "\t" + glycans.Count + "\t" + glycansSeenOnThisProtein);
            }


            ////////////////////////////////////////////////////////////count and print glycopeptides output

            outputPeptides.WriteLine("PeptideWithMods\tPeptide\tMods\tuniprotID\tProteinName\tLocalized\t#ofPSMs\t#ofLocalizedPSMs\t#OfGlycoSitesInSeq\t#Glycans" +
                    "\tGlycansSeen\tGlycanTypes\t#LocalizedGlycans\tLocalizedGlycansSeen\tLocalizedGlycoTypes\tLinkage" +
                    "\tNXS\tNXT\t#ofChargeStates\tChargeStatesSeen\tAvgScore\tStdDevScore\tMedianScore\tBestScore\tfragOfBestScore\tseenInUniprot\tevidenceType" +
                    "\tevidenceNumber\tSeenWithHCD\tSeenWithETD\tSeenWithBoth\tAvgGlcNAcRatio\tStdDevGlcNAcRatio\tmedianGlcNAcRatio\tAvgRT\tStdDevRT" +
                    "\tDeltaRTClosestPSM\tDeltaRTFurthestPSM\tAvgDeltaRT\tStdDevDeltaRT\tCharge\tPrecursor Theoretical m/z (Th)\tSpectrum number\tLFQ Intensity\tByonic Intensity");

            int numberOfUniqueGlycopeptides = glycoPeptideDictionary.Count;
            int numberOfUniqueGlycopeptidesLocalized = 0;
            int numberOfUniqueGlycopeptidesNlinked = 0;
            int numberOfUniqueGlycopeptidesNlinkedLocalized = 0;
            int numberOfUniqueGlycopeptidesOlinked = 0;
            int numberOfUniqueGlycopeptidesOlinkedLocalized = 0;

            int numberOfUniqueETD = 0;
            int numberOfUniqueETDLocalized = 0;
            int numberOfUniqueETD_Nlinked = 0;
            int numberOfUniqueETD_NlinkedLocalized = 0;
            int numberOfUniqueETD_Olinked = 0;
            int numberOfUniqueETD_OlinkedLocalized = 0;
            int numberOfUniqueHCD = 0;
            int numberOfHCDLocalized = 0;
            int numberOfUniqueHCD_Nlinked = 0;
            int numberOfUniqueHCD_NlinkedLocalized = 0;
            int numberOfUniqueHCD_Olinked = 0;
            int numberOfUniqueHCD_OlinkedLocalized = 0;

            double numberOfTimesSeenWithBoth = 0;
            double numberOfTimesETDbetter = 0;

            foreach (KeyValuePair<string, List<GlycoPSM>> entry in glycoPeptideDictionary)
            {
                bool Nlinked = false;
                bool NlinkedLocalized = false;
                bool Olinked = false;
                bool OlinkedLocalized = false;

                bool etdLocalized = false;
                bool etdNlinked = false;
                bool etdNlinkedLocalized = false;
                bool etdOlinked = false;
                bool etdOlinkedLocalized = false;
                bool hcdLocalized = false;
                bool hcdNlinked = false;
                bool hcdNlinkedLocalized = false;
                bool hcdOlinked = false;
                bool hcdOlinkedLocalized = false;

                int numberOfPSMs = 0;
                int numberOfPSMsLocalized = 0;

                string sequenceToUse = entry.Key.Split(' ')[0];

                List<string> mods = new List<string>();

                string uniprotID = "none";
                string proteinName = "none";
                string peptideSequenceOnly = "none";

                bool isLocalized = false;
                int numberOfGlycoSites = 0;

                HashSet<string> glycansInThisPSM = new HashSet<string>();
                HashSet<string> glycanTypesInThisPSM = new HashSet<string>();
                HashSet<string> glycansInThisPSMlocalized = new HashSet<string>();
                HashSet<string> glycanTypesInThisPSMLocalized = new HashSet<string>();


                bool NXS = false;
                bool NXT = false;

                string glycotype = "unAssigned";

                List<double> retentionTimes = new List<double>();

                HashSet<int> chargeState = new HashSet<int>();
                List<double> scores = new List<double>();

                double bestScore = 0;
                string fragmentationOfBestScore = "notAssigned";
                GlycoPSM bestScoringPSM = entry.Value[0];

                bool seenInUniprot = false;
                string evidenceType = "None";
                int evidenceNumber = 0;

                bool seenWithHCD = false;
                bool seenWithETD = false;
                bool seenWithBoth = false;

                List<double> glcNAcToGalNAcRatios = new List<double>();

                foreach (GlycoPSM psm in entry.Value)
                {
                    uniprotID = psm.uniprotID;
                    proteinName = psm.proteinName;
                    peptideSequenceOnly = psm.peptide.Sequence;
                    numberOfGlycoSites = psm.numberOfSites;

                    mods = psm.mods;

                    numberOfPSMs++;
                    foreach (string glycan in psm.glycans)
                    {
                        glycansInThisPSM.Add(glycan);
                    }
                    
                    foreach(string glycanType in psm.glycanTypes)
                    {
                        glycanTypesInThisPSM.Add(glycanType);
                    }

                    if (psm.Nlinked && psm.Olinked)
                        glycotype = "BothNandO";
                    else if (psm.Nlinked)
                    {
                        glycotype = "Nlinked";
                        if (psm.NXSmotif)
                            NXS = true;
                        if (psm.NXTmotif)
                            NXT = true;
                        Nlinked = true;
                        if (psm.seenWithETD)
                            etdNlinked = true;
                        if (psm.seenWithHCD)
                            hcdNlinked = true;
                    }
                    else if (psm.Olinked)
                    {
                        glycotype = "Olinked";
                        Olinked = true;
                        if (psm.seenWithETD)
                            etdOlinked = true;
                        if (psm.seenWithHCD)
                            hcdOlinked = true;
                    }

                    chargeState.Add(psm.charge);
                    scores.Add(psm.score);

                    // Get best scoring psm
                    if (psm.score > bestScore)
                    {
                        bestScore = psm.score;
                        fragmentationOfBestScore = psm.fragmentation;
                        bestScoringPSM = psm;
                    }

                    if (psm.matchedToUniprot)
                    {
                        seenInUniprot = true;
                        evidenceNumber = psm.evidenceNumber;
                        evidenceType = psm.evidenceType;
                    }

                    if (psm.fragmentation.Equals("HCD"))
                        seenWithHCD = true;
                    else if (psm.fragmentation.Equals("ETD"))
                        seenWithETD = true;

                    retentionTimes.Add(psm.scanTime);
                    glcNAcToGalNAcRatios.Add(psm.GlcNAcGalNAcratio);

                    if (psm.isLocalized)
                    {
                        isLocalized = true;
                        numberOfPSMsLocalized++;

                        foreach (string glycan in psm.glycans)
                        {
                            glycansInThisPSMlocalized.Add(glycan);
                        }

                        foreach(string glycanType in psm.glycanTypes)
                        {
                            glycanTypesInThisPSMLocalized.Add(glycanType);
                        }

                        if (psm.Nlinked)
                        {
                            NlinkedLocalized = true;
                            if (psm.seenWithETD)
                                etdNlinkedLocalized = true;
                            if (psm.seenWithHCD)
                                hcdNlinkedLocalized = true;
                        }
                        else if (psm.Olinked)
                        {
                            OlinkedLocalized = true;
                            if (psm.seenWithETD)
                                etdOlinkedLocalized = true;
                            if (psm.seenWithHCD)
                                hcdOlinkedLocalized = true;
                        }

                        if (psm.fragmentation.Equals("HCD"))
                            hcdLocalized = true;
                        else if (psm.fragmentation.Equals("ETD"))
                            etdLocalized = true;

                    }
                }

                double avgScore = Statistics.Mean(scores);
                double stdDevScore = Statistics.StandardDeviation(scores);
                double medianScore = Statistics.Median(scores);

                double avgGlcNAcRatio = Statistics.Mean(glcNAcToGalNAcRatios);
                double stdDevGlcNAcRatio = Statistics.StandardDeviation(glcNAcToGalNAcRatios);
                double medianGlcNAcRatio = Statistics.Median(glcNAcToGalNAcRatios);

                if (seenWithHCD && seenWithETD)
                {
                    seenWithBoth = true;
                    numberOfTimesSeenWithBoth++;
                }

                if (fragmentationOfBestScore.Equals("ETD") && seenWithBoth)
                {
                    numberOfTimesETDbetter++;
                }

                string glycans = "none";
                foreach (string glycanToAdd in glycansInThisPSM)
                {
                    if (glycans.Equals("none"))
                        glycans = glycanToAdd;
                    else
                        glycans += ";" + glycanToAdd;
                }

                string glycanTypes = "none";
                foreach(string glycanType in glycanTypesInThisPSM)
                {
                    if (glycanTypes.Equals("none"))
                        glycanTypes = glycanType;
                    else
                        glycanTypes += ";" + glycanType;
                }


                string glycansLocalized = "none";
                foreach (string glycanToAdd in glycansInThisPSMlocalized)
                {
                    if (glycansLocalized.Equals("none"))
                        glycansLocalized = glycanToAdd;
                    else
                        glycansLocalized += ";" + glycanToAdd;
                }

                string glycanTypesLocalized = "none";
                foreach (string glycanType in glycanTypesInThisPSM)
                {
                    if (glycanTypesLocalized.Equals("none"))
                        glycanTypesLocalized = glycanType;
                    else
                        glycanTypesLocalized += ";" + glycanType;
                }

                int numberOfChargeStates = chargeState.Count;
                string chargeStatesSeen = "none";
                foreach (int z in chargeState)
                {
                    if (chargeStatesSeen.Equals("none"))
                        chargeStatesSeen = z.ToString();
                    else
                        chargeStatesSeen += ";" + z.ToString();
                }

                string modsOnPep = "none";
                foreach (string mod in mods)
                {
                    if (modsOnPep.Equals("none"))
                        modsOnPep = mod;
                    else
                        modsOnPep += ";" + mod;
                }

                double avgRT = Statistics.Mean(retentionTimes);
                double stdDevRT = Statistics.StandardDeviation(retentionTimes);

                double deltaRTbetweenClosestPSM = GetDeltaOfClosestRT(retentionTimes);
                double deltaRTbetweenFurthestPSM = GetDeltaOfFurthestRT(retentionTimes);
                double avgDeltaRT = GetAvgDeltaRT(retentionTimes);
                double stdDevDeltaRT = GetStdDevDeltaRT(retentionTimes);

                if (isLocalized)
                    numberOfUniqueGlycopeptidesLocalized++;
                if (Nlinked)
                    numberOfUniqueGlycopeptidesNlinked++;
                if (NlinkedLocalized)
                    numberOfUniqueGlycopeptidesNlinkedLocalized++;
                if (Olinked)
                    numberOfUniqueGlycopeptidesOlinked++;
                if (OlinkedLocalized)
                    numberOfUniqueGlycopeptidesOlinkedLocalized++;

                if (seenWithHCD)
                    numberOfUniqueHCD++;
                if (hcdLocalized)
                    numberOfHCDLocalized++;
                if (hcdNlinked)
                    numberOfUniqueHCD_Nlinked++;
                if (hcdNlinkedLocalized)
                    numberOfUniqueHCD_NlinkedLocalized++;
                if (hcdOlinked)
                    numberOfUniqueHCD_Olinked++;
                if (hcdOlinkedLocalized)
                    numberOfUniqueHCD_OlinkedLocalized++;

                if (seenWithETD)
                    numberOfUniqueETD++;
                if (etdLocalized)
                    numberOfUniqueETDLocalized++;
                if (etdNlinked)
                    numberOfUniqueETD_Nlinked++;
                if (etdNlinkedLocalized)
                    numberOfUniqueETD_NlinkedLocalized++;
                if (etdOlinked)
                    numberOfUniqueETD_Olinked++;
                if (etdOlinkedLocalized)
                    numberOfUniqueETD_OlinkedLocalized++;
                

                outputPeptides.WriteLine(sequenceToUse + "\t" + peptideSequenceOnly + "\t" + modsOnPep + "\t" + uniprotID + "\t" + proteinName + "\t" + isLocalized + "\t" + numberOfPSMs + "\t" + numberOfPSMsLocalized + "\t" + numberOfGlycoSites + "\t" +
                    glycansInThisPSM.Count + "\t" + glycans + "\t" + glycanTypes + "\t" + glycansInThisPSMlocalized.Count + "\t" + glycansLocalized + "\t" + glycanTypesLocalized + "\t" + glycotype + "\t" + NXS + "\t" + NXT + "\t" + numberOfChargeStates + "\t" + chargeStatesSeen + "\t" + avgScore + "\t" +
                    stdDevScore + "\t" + medianScore + "\t" + bestScore + "\t" + fragmentationOfBestScore + "\t" + seenInUniprot + "\t" + evidenceType + "\t" + evidenceNumber + "\t" + seenWithHCD + "\t" + seenWithETD + "\t" +
                    seenWithBoth + "\t" + avgGlcNAcRatio + "\t" + stdDevGlcNAcRatio + "\t" + medianGlcNAcRatio + "\t" + avgRT + "\t" + stdDevRT + "\t" + deltaRTbetweenClosestPSM + "\t" + deltaRTbetweenFurthestPSM + "\t" + avgDeltaRT + "\t" + stdDevDeltaRT +
                    "\t" + bestScoringPSM.charge + "\t" + bestScoringPSM.mzCalc + "\t" + bestScoringPSM.scanNumber + "\t" + bestScoringPSM.intensity + "\t" + bestScoringPSM.byonicIntensity);
            }


            ////////////////////////////////////////////////////////////count and print uniqueSequences

            outputSequences.WriteLine("Sequence\tMods\tuniprotID\tProteinName\tLocalized\t#ofPSMs\t#ofLocalizedPSMs\t#UniqPeps\t#UniqPepsLocalized\t#UniqGlycoIsoforms\t#UniqGlycoIsoformsLocalized\t#OfGlycoSitesInSeq\t#Glycans\tGlycansSeen" +
                    "\t#LocalizedGlycans\tLocalizedGlycansSeen\tGlycoType\tNXS\tNXT\t#ofChargeStates\tChargeStatesSeen\tPepMassNoMods\tAvgScore\tStdDevScore\tMedianScore\tBestScore\tfragOfBestScore\tseenInUniprot\tevidenceType\tevidenceNumber\tSeenWithHCD" +
                    "\tSeenWithETD\tSeenWithBoth\tAvgGlcNAcRatio\tStdDevGlcNAcRatio\tmedianGlcNAcRatio\tAvgRT\tStdDevRT\tDeltaRTClosestPSM\tDeltaRTFurthestPSM\tAvgDeltaRT\tStdDevDeltaRT");

            int countOfUniqueSequences = uniqueSequencesDictionary.Count;
            int countOfUniqueSequencesLocalized = 0;
            int countOfUniqueSequencesNlinked = 0;
            int countOfUniqueSequencesNlinkedLocalized = 0;
            int countOfUniqueSequencesOlinked = 0;
            int countOfUniqueSequencesOlinkedLocalized = 0;

            int uniqueSeqWithNXS = 0;
            int uniqueSeqWithNXSlocalized = 0;
            int uniqueSeqWithNXT = 0;
            int uniqueSeqWithNXTlocalized = 0;


            foreach (KeyValuePair<string, List<GlycoPSM>> entry in uniqueSequencesDictionary)
            {
                bool Nlinked = false;
                bool NlinkedLocalized = false;
                bool Olinked = false;
                bool OlinkedLocalized = false;

                int numberOfPSMs = 0;
                int numberOfPSMsLocalized = 0;

                HashSet<string> uniquePeps = new HashSet<string>();
                HashSet<string> uniquePepsLocalized = new HashSet<string>();
                HashSet<string> uniqueGlycoIsoforms = new HashSet<string>();
                HashSet<string> uniqueGlycoIsoformsLocalized = new HashSet<string>();


                string sequenceToUse = entry.Key;

                List<string> mods = new List<string>();

                string uniprotID = "none";
                string proteinName = "none";

                bool isLocalized = false;
                int numberOfGlycoSites = 0;

                HashSet<string> glycansInThisSequence = new HashSet<string>();
                HashSet<string> glycansInThisSequenceLocalized = new HashSet<string>();

                bool NXS = false;
                bool NXT = false;

                bool NXSlocalized = false;
                bool NXTlocalized = false;

                string glycotype = "unAssigned";

                List<double> retentionTimes = new List<double>();

                HashSet<int> chargeState = new HashSet<int>();
                List<double> scores = new List<double>();

                double bestScore = 0;
                string fragmentationOfBestScore = "notAssigned";

                bool seenInUniprot = false;
                string evidenceType = "None";
                int evidenceNumber = 0;

                bool seenWithHCD = false;
                bool seenWithETD = false;
                bool seenWithBoth = false;

                double peptideMassNoMods = 0;

                List<double> glcNAcToGalNAcRatios = new List<double>();

                foreach (GlycoPSM psm in entry.Value)
                {
                    uniprotID = psm.uniprotID;
                    proteinName = psm.proteinName;
                    numberOfGlycoSites = psm.numberOfSites;

                    peptideMassNoMods = psm.peptideMonoMass;

                    uniquePeps.Add(psm.sequenceWithMods);

                    mods = psm.mods;

                    numberOfPSMs++;
                    string glycansForIsoform = "none";

                    foreach (string glycan in psm.glycans)
                    {
                        glycansInThisSequence.Add(glycan);
                        if (glycansForIsoform.Equals("none"))
                            glycansForIsoform = glycan;
                        else
                            glycansForIsoform += "_" + glycan;
                    }

                    string glycoIsoform = psm.peptide.Sequence + glycansForIsoform;
                    uniqueGlycoIsoforms.Add(glycoIsoform);

                    if (psm.Nlinked && psm.Olinked)
                        glycotype = "BothNandO";
                    else if (psm.Nlinked)
                    {
                        glycotype = "Nlinked";
                        if (psm.NXSmotif)
                            NXS = true;
                        if (psm.NXTmotif)
                            NXT = true;
                        Nlinked = true;
                    }
                    else if (psm.Olinked)
                    {
                        glycotype = "Olinked";
                        Olinked = true;
                    }

                    chargeState.Add(psm.charge);
                    scores.Add(psm.score);

                    if (psm.score > bestScore)
                    {
                        bestScore = psm.score;
                        fragmentationOfBestScore = psm.fragmentation;
                    }

                    if (psm.matchedToUniprot)
                    {
                        seenInUniprot = true;
                        evidenceNumber = psm.evidenceNumber;
                        evidenceType = psm.evidenceType;
                    }

                    if (psm.fragmentation.Equals("HCD"))
                        seenWithHCD = true;
                    else if (psm.fragmentation.Equals("ETD"))
                        seenWithETD = true;

                    retentionTimes.Add(psm.scanTime);
                    glcNAcToGalNAcRatios.Add(psm.GlcNAcGalNAcratio);

                    if (psm.isLocalized)
                    {
                        isLocalized = true;
                        numberOfPSMsLocalized++;

                        uniquePepsLocalized.Add(psm.sequenceWithMods);

                        string glycansForIsoformLocalized = "none";

                        foreach (string glycan in psm.glycans)
                        {
                            glycansInThisSequenceLocalized.Add(glycan);
                            if (glycansForIsoformLocalized.Equals("none"))
                                glycansForIsoformLocalized = glycan;
                            else
                                glycansForIsoformLocalized += "_" + glycan;
                        }

                        string glycoIsoformLocalized = psm.peptide.Sequence + glycansForIsoformLocalized;
                        uniqueGlycoIsoformsLocalized.Add(glycoIsoformLocalized);

                        if (psm.Nlinked)
                        {
                            NlinkedLocalized = true;
                        }
                        else if (psm.Olinked)
                        {
                            OlinkedLocalized = true;
                        }

                        if (psm.NXSmotif)
                            NXSlocalized = true;

                        if (psm.NXTmotif)
                            NXTlocalized = true;

                    }


                }

                double avgScore = Statistics.Mean(scores);
                double stdDevScore = Statistics.StandardDeviation(scores);
                double medianScore = Statistics.Median(scores);

                double avgGlcNAcRatio = Statistics.Mean(glcNAcToGalNAcRatios);
                double stdDevGlcNAcRatio = Statistics.StandardDeviation(glcNAcToGalNAcRatios);
                double medianGlcNAcRatio = Statistics.Median(glcNAcToGalNAcRatios);

                if (seenWithHCD && seenWithETD)
                {
                    seenWithBoth = true;
                    numberOfTimesSeenWithBoth++;
                }

                if (fragmentationOfBestScore.Equals("ETD") && seenWithBoth)
                {
                    numberOfTimesETDbetter++;
                }

                string glycans = "none";
                foreach (string glycanToAdd in glycansInThisSequence)
                {
                    if (glycans.Equals("none"))
                        glycans = glycanToAdd;
                    else
                        glycans += ";" + glycanToAdd;
                }


                string glycansLocalized = "none";
                foreach (string glycanToAdd in glycansInThisSequenceLocalized)
                {
                    if (glycansLocalized.Equals("none"))
                        glycansLocalized = glycanToAdd;
                    else
                        glycansLocalized += ";" + glycanToAdd;
                }

                int numberOfChargeStates = chargeState.Count;
                string chargeStatesSeen = "none";
                foreach (int z in chargeState)
                {
                    if (chargeStatesSeen.Equals("none"))
                        chargeStatesSeen = z.ToString();
                    else
                        chargeStatesSeen += ";" + z.ToString();
                }

                string modsOnPep = "none";
                foreach (string mod in mods)
                {
                    if (modsOnPep.Equals("none"))
                        modsOnPep = mod;
                    else
                        modsOnPep += ";" + mod;
                }

                double avgRT = Statistics.Mean(retentionTimes);
                double stdDevRT = Statistics.StandardDeviation(retentionTimes);

                double deltaRTbetweenClosestPSM = GetDeltaOfClosestRT(retentionTimes);
                double deltaRTbetweenFurthestPSM = GetDeltaOfFurthestRT(retentionTimes);
                double avgDeltaRT = GetAvgDeltaRT(retentionTimes);
                double stdDevDeltaRT = GetStdDevDeltaRT(retentionTimes);

                if (isLocalized)
                    countOfUniqueSequencesLocalized++;
                if (Nlinked)
                    countOfUniqueSequencesNlinked++;
                if (NlinkedLocalized)
                    countOfUniqueSequencesNlinkedLocalized++;
                if (Olinked)
                    countOfUniqueSequencesOlinked++;
                if (OlinkedLocalized)
                    countOfUniqueSequencesOlinkedLocalized++;

                if (NXS)
                    uniqueSeqWithNXS++;
                if (NXT)
                    uniqueSeqWithNXT++;
                if (NXSlocalized)
                    uniqueSeqWithNXSlocalized++;
                if (NXTlocalized)
                    uniqueSeqWithNXTlocalized++;

                outputSequences.WriteLine(sequenceToUse + "\t" + modsOnPep + "\t" + uniprotID + "\t" + proteinName + "\t" + isLocalized + "\t" + numberOfPSMs + "\t" + numberOfPSMsLocalized + "\t" + uniquePeps.Count + "\t" + uniquePepsLocalized.Count + "\t" +
                    uniqueGlycoIsoforms.Count + "\t" + uniqueGlycoIsoformsLocalized.Count + "\t" + numberOfGlycoSites + "\t" + glycansInThisSequence.Count + "\t" + glycans + "\t" + glycansInThisSequenceLocalized.Count + "\t" + glycansLocalized + "\t" +
                    glycotype + "\t" + NXS + "\t" + NXT + "\t" + numberOfChargeStates + "\t" + chargeStatesSeen + "\t" + peptideMassNoMods + "\t" + avgScore + "\t" + stdDevScore + "\t" + medianScore + "\t" + bestScore + "\t" + fragmentationOfBestScore + "\t" + seenInUniprot + "\t" +
                    evidenceType + "\t" + evidenceNumber + "\t" + seenWithHCD + "\t" + seenWithETD + "\t" + seenWithBoth + "\t" + avgGlcNAcRatio + "\t" + stdDevGlcNAcRatio + "\t" + medianGlcNAcRatio + "\t" + avgRT + "\t" + stdDevRT + "\t" +
                    deltaRTbetweenClosestPSM + "\t" + deltaRTbetweenFurthestPSM + "\t" + avgDeltaRT + "\t" + stdDevDeltaRT);
            }


            ////////////////////////////////////////////////////////////count and print PSMs

            outputPSMs.WriteLine("Sequence\tUniprot\tProteinName\tseqWithMods\tPeptideStartPos.\tMods\tGlycans\tGlycanTypes\tGlycanPositions\t#GlycoSitesOnPep\t"+
                                "PEP2D\tFDR2D\tFDR2DUniq\tqvalue2D\tPEP1D\tFDR1D\tFDR1DUniq\tqvalue1D\tlogProb\tScore\tDeltaScore\tDeltaModScore\tCharge\tmzObs"+
                                "\tmxCalc\tobsMH\tcalcMH\tPepMassNoMod\tCleavage\tscanTime\tScanNumber\tMasterScan\tFragmentation\tIsGlycoPeptide\tSeenWithHCD\t"+
                                "SeenWithETD\tLocalized\tNlinked\tOlinked\tNXSmotif\tNXTmotif\tMatchedToUniprot\tEvidence\tEvidence#\tGlcNAcRatio\tFile");

            int countPSMs = 0;
            int countLocalizedPSMs = 0;
            int countNlinkedPSMs = 0;
            int countNlinkedPSMsLocalized = 0;
            int countOlinkedPSMs = 0;
            int countOlinkedPSMsLocalized = 0;

            int psmsETD = 0;
            int psmsETDLocalized = 0;
            int psmsETD_Nlinked = 0;
            int psmsETD_NlinkedLocalized = 0;
            int psmsETD_Olinked = 0;
            int psmsETD_OlinkedLocalized = 0;
            int psmsHCD = 0;
            int psmsHCDLocalized = 0;
            int psmsHCD_Nlinked = 0;
            int psmsHCD_NlinkedLocalized = 0;
            int psmsHCD_Olinked = 0;
            int psmsHCD_OlinkedLocalized = 0;

            foreach (GlycoPSM psm in glycoPSMs)
            {
                if (!psm.proteinName.Contains(">Reverse") && !psm.proteinName.Contains("DECOY"))
                {
                    countPSMs++;
                    if (psm.seenWithHCD)
                    {
                        psmsHCD++;
                        if (psm.Nlinked)
                        {
                            psmsHCD_Nlinked++;
                        }
                        if (psm.Olinked)
                        {
                            psmsHCD_Olinked++;
                        }
                    }
                    if (psm.seenWithETD)
                    {
                        psmsETD++;
                        if (psm.Nlinked)
                        {
                            psmsETD_Nlinked++;
                        }
                        if (psm.Olinked)
                        {
                            psmsETD_Olinked++;
                        }
                    }

                    if (psm.Nlinked)
                        countNlinkedPSMs++;
                    if (psm.Olinked)
                        countOlinkedPSMs++;

                    if (psm.isLocalized)
                    {
                        countLocalizedPSMs++;
                        if (psm.seenWithHCD)
                        {
                            psmsHCDLocalized++;
                            if (psm.Nlinked)
                            {
                                psmsHCD_NlinkedLocalized++;
                            }
                            if (psm.Olinked)
                            {
                                psmsHCD_OlinkedLocalized++;
                            }
                        }
                        if (psm.seenWithETD)
                        {
                            psmsETDLocalized++;
                            if (psm.Nlinked)
                            {
                                psmsETD_NlinkedLocalized++;
                            }
                            if (psm.Olinked)
                            {
                                psmsETD_OlinkedLocalized++;
                            }
                        }

                        if (psm.Nlinked)
                            countNlinkedPSMsLocalized++;
                        if (psm.Olinked)
                            countOlinkedPSMsLocalized++;
                    }

                    /*
                     *  GlycoPSM psm = new GlycoPSM(peptide, peptide.SequenceWithModifications, mods, glycans, glycanPositionsList, uniprotID, PEP2D, PEP1D, logProb, score, deltaScore, deltaModScore, mzObs, mzCalc, charge, numberOfSites, ppmError,
                            obsMH, calcMH, cleavage, proteinName, peptideStartPosition, scanTime, scanNumber, FDR2D, FDR1D, FDR2Dunique, FDR1Dunique, qvalue2D, qvalue1D, fragmentation, masterScan, peak126, peak138, peak144,
                            peak168, peak186, peak204, peak274, peak292, peak366, GlcNAcGalNAcratio, isGlycoPeptide, seenWithHCD, seenWithETD, NXSmotif, NXTmotif, isLocalized, Nlinked, Olinked, matchedToUniprot, uniprotEvidenceType, uniprotEvidenceNumber);
                     * */

                    string modsOnPep = "none";
                    foreach (string mod in psm.mods)
                    {
                        if (modsOnPep.Equals("none"))
                            modsOnPep = mod;
                        else
                            modsOnPep += ";" + mod;
                    }

                    string glycanTypes = "none";
                    foreach(string type in psm.glycanTypes)
                    {

                        if (glycanTypes.Equals("none"))
                            glycanTypes = type;
                        else
                            glycanTypes += ";" + type;
                    }

                    string glycansOnPep = "none";
                    foreach (string glycan in psm.glycans)
                    {
                        if (glycansOnPep.Equals("none"))
                            glycansOnPep = glycan;
                        else
                            glycansOnPep += ";" + glycan;
                    }

                    string glycanPositionsInPep = "none";
                    foreach (int position in psm.glycanPositions)
                    {
                        if (glycanPositionsInPep.Equals("none"))
                            glycanPositionsInPep = position.ToString();
                        else
                            glycanPositionsInPep += ";" + position;
                    }



                    outputPSMs.WriteLine(psm.peptide.Sequence + "\t" + psm.uniprotID + "\t" + psm.proteinName + "\t" + psm.sequenceWithMods + "\t" + 
                        psm.peptideStartPosition + "\t" + modsOnPep + "\t" + glycansOnPep + "\t" + glycanTypes + "\t" + glycanPositionsInPep + "\t" + 
                        psm.numberOfSites + "\t" + psm.PEP2D + "\t" + psm.FDR2D + "\t" + psm.FDR2Dunique + "\t" + psm.qvalue2D + "\t" + psm.PEP1D + 
                        "\t" + psm.FDR1D + "\t" + psm.FDR1Dunique + "\t" + psm.qvalue1D + "\t" + psm.logProb + "\t" + psm.score + "\t" + psm.deltaScore + 
                        "\t" + psm.deltaModScore  + "\t" + psm.charge + "\t" + psm.mzObs + "\t" + psm.mzCalc + "\t" + psm.obsMH + "\t" + psm.calcMH + "\t" + 
                        psm.peptideMonoMass + "\t" + psm.cleavage + "\t" + psm.scanTime + "\t" + psm.scanNumber + "\t" + psm.masterScan + "\t" + 
                        psm.fragmentation + "\t" + psm.isGlycoPeptide + "\t" + psm.seenWithHCD + "\t" + psm.seenWithETD + "\t" + psm.isLocalized + "\t" 
                        + psm.Nlinked + "\t" + psm.Olinked + "\t" + psm.NXSmotif + "\t" + psm.NXTmotif + "\t" + psm.matchedToUniprot + "\t" + psm.evidenceType + 
                        "\t" + psm.evidenceNumber + "\t" + psm.GlcNAcGalNAcratio + "\t" + psm.file);
                }
            }

            OnUpdateProgress(1.0);

            ////////////////////////////////////////////////////////////write summary

            outputSummary.WriteLine("File:\t" + outputPath);
            outputSummary.WriteLine("Localized\tTotal\tNlinked\tOlinked");
            outputSummary.WriteLine("Localized GlycoPSMs" + "\t" + countLocalizedPSMs + "\t" + countNlinkedPSMsLocalized + "\t" + countOlinkedPSMsLocalized);
            outputSummary.WriteLine("Localized HCD PSMs" + "\t" + psmsHCDLocalized + "\t" + psmsHCD_NlinkedLocalized + "\t" + psmsHCD_OlinkedLocalized);
            outputSummary.WriteLine("Localized ETD PSMs" + "\t" + psmsETDLocalized + "\t" + psmsETD_NlinkedLocalized + "\t" + psmsETD_OlinkedLocalized);
            outputSummary.WriteLine("Localized GlycoPeps" + "\t" + numberOfUniqueGlycopeptidesLocalized + "\t" + numberOfUniqueGlycopeptidesNlinkedLocalized + "\t" + numberOfUniqueGlycopeptidesOlinkedLocalized);
            outputSummary.WriteLine("Localized HCD Peps" + "\t" + numberOfHCDLocalized + "\t" + numberOfUniqueHCD_NlinkedLocalized + "\t" + numberOfUniqueHCD_OlinkedLocalized);
            outputSummary.WriteLine("Localized ETD Peps" + "\t" + numberOfUniqueETDLocalized + "\t" + numberOfUniqueETD_NlinkedLocalized + "\t" + numberOfUniqueETD_OlinkedLocalized);
            outputSummary.WriteLine("Localized UniqueSeq" + "\t" + countOfUniqueSequencesLocalized + "\t" + countOfUniqueSequencesNlinkedLocalized + "\t" + countOfUniqueSequencesOlinkedLocalized);
            outputSummary.WriteLine("Localized GlycoProteins" + "\t" + countOfGlycoProteinsLocalized + "\t" + countOfGlycoProteinsNlinkedLocalized + "\t" + countOfGlycoProteinsOlinkedLocalized);
            outputSummary.WriteLine("Localized GlycoSites" + "\t" + countOfGlycoSitesLocalized + "\t" + countOfGlycoSitesLocazliedNlinked + "\t" + countOfGlycoSitesLocalizedOlinked);
            outputSummary.WriteLine("Localized Glycans" + "\t" + countOfGlycansLocalized + "\t" + countOfGlycansLocazliedNlinked + "\t" + countOfGlycansLocalizedOlinked);
            outputSummary.WriteLine("All\tTotal\tNlinked\tOlinked");
            outputSummary.WriteLine("All GlycoPSMs" + "\t" + countPSMs + "\t" + countNlinkedPSMs + "\t" + countOlinkedPSMs);
            outputSummary.WriteLine("All HCD PSMs" + "\t" + psmsHCD + "\t" + psmsHCD_Nlinked + "\t" + psmsHCD_Olinked);
            outputSummary.WriteLine("All ETD PSMs" + "\t" + psmsETD + "\t" + psmsETD_Nlinked + "\t" + psmsETD_Olinked);
            outputSummary.WriteLine("All GlycoPeps" + "\t" + numberOfUniqueGlycopeptides + "\t" + numberOfUniqueGlycopeptidesNlinked + "\t" + numberOfUniqueGlycopeptidesOlinked);
            outputSummary.WriteLine("All HCD Peps" + "\t" + numberOfUniqueHCD + "\t" + numberOfUniqueHCD_Nlinked + "\t" + numberOfUniqueHCD_Olinked);
            outputSummary.WriteLine("All ETD Peps" + "\t" + numberOfUniqueETD + "\t" + numberOfUniqueETD_Nlinked + "\t" + numberOfUniqueETD_Olinked);
            outputSummary.WriteLine("All UniqueSeq" + "\t" + countOfUniqueSequences + "\t" + countOfUniqueSequencesNlinked + "\t" + countOfUniqueSequencesOlinked);
            outputSummary.WriteLine("All GlycoProteins" + "\t" + countOfGlycoProteins + "\t" + countOfGlycoProteinsNlinked + "\t" + countOfGlycoProteinsOlinked);
            outputSummary.WriteLine("All GlycoSites" + "\t" + countOfGlycoSitesAll + "\t" + countOfGlycoSitesAllNlinked + "\t" + countOfGlycoSitesAllOlinked);
            outputSummary.WriteLine("All Glycans" + "\t" + countOfGlycansAll + "\t" + countOfGlcyansAllNlinked + "\t" + countOfGlycansAllOlinked);


            /**
            outputSummary.WriteLine();
            outputSummary.WriteLine("\tSeenWithBothHCD+ETD\tETDisBetterScore\tPercentage");
            outputSummary.WriteLine("UniqueGlycoPeps\t" + numberOfTimesSeenWithBoth + "\t" + numberOfTimesETDbetter + "\t" + (numberOfTimesETDbetter / numberOfTimesSeenWithBoth) * 100);
            **/

            outputSummary.WriteLine();
            outputSummary.WriteLine("\tOverall FDR\tHCD FDR\tETD FDR");
            outputSummary.WriteLine("FDR at Localized GlycoPSM level Filtered\t" + fdrAllLocalized + "\t" + fdrHCDLocalized + "\t" + fdrETDLocalized);
            outputSummary.WriteLine("FDR at GlycoPSM level Filtered\t" + fdrAllGlycoPSMs + "\t" + fdrHCD + "\t" + fdrETD);
            outputSummary.WriteLine("FDR at GlycoPSM level Raw\t" + fdrAllGlycoPSMsNoFilter + "\t" + fdrHCDNoFilter + "\t" + fdrETDNoFilter);

            outputSummary.WriteLine();
            outputSummary.WriteLine("Sites In Uniprot\tLocalized\tAll");
            outputSummary.WriteLine("Total\t" + countOfSitesSeenInUniprotLocalized + "\t" + countOfSitesSeenInUniprot);
            outputSummary.WriteLine("Report in Lit\t" + countReportInLitLocalized + "\t" + countReportInLit);
            outputSummary.WriteLine("Sequence Analysis\t" + countSeqAnalysisLocalized + "\t" + countSeqAnalysis);
            outputSummary.WriteLine("By Similarity\t" + countBySimilarityLocalized + "\t" + countBySimilarity);
            outputSummary.WriteLine("Inferred By Curator\t" + countInferredByCuratorLocalized + "\t" + countInferredByCurator);
            outputSummary.WriteLine("Inferred By Evidence\t" + countInferredEvidenceLocalized + "\t" + countInferredEvidence);
            outputSummary.WriteLine("Discussed, No Evidence\t" + countDiscussedNoEvidenceLocalized + "\t" + countDiscussedNoEvidence);

            outputSummary.WriteLine();
            outputSummary.WriteLine("Motifs Analysis (Uniq.Seq.)\tLocalized\tAll");
            outputSummary.WriteLine("N-X-S\t" + uniqueSeqWithNXSlocalized + "\t" + uniqueSeqWithNXS);
            outputSummary.WriteLine("N-X-T\t" + uniqueSeqWithNXTlocalized + "\t" + uniqueSeqWithNXT);

            outputSummary.WriteLine();
            outputSummary.WriteLine("Notes:");
            outputSummary.WriteLine("1.)\tAll results here are filtered to 1% FDR at the PSM level using Byonic's 2D-FDR value");
            outputSummary.WriteLine("\tand glycopeptides must be at least 5 amino acids long to be considered.");
            outputSummary.WriteLine("\tThe 'Raw' FDR calculation above shows the FDR if these two filters are not applied.");
            outputSummary.WriteLine("2.)\tThe Uniprot Database used to match Glycosites was downloaded April 10, 2017 as a gff file");
            outputSummary.WriteLine("\tand was compiled using the 'UniProtGlycoSiteParser' program written by NMR.");
            outputSummary.WriteLine("3.)\tThere are seven files printed through the GlycoCompiler2 program, including this summary");
            outputSummary.WriteLine("\tand files for GlycoPSMs, GlycoPeptides, UniqueSequences, GlycoProteins, GlycoSites, and Glycans.");


            printStats(glycoPSMs, @outputFolder + "\\_stats.csv", glycanMasses);



            ////////////////////////////////////////////////////////////close streamwriters
            outputSites.Close();
            outputGlycans.Close();
            outputProteins.Close();
            outputSummary.Close();
            outputPeptides.Close();
            outputSequences.Close();
            outputPSMs.Close();
            OnUpdateProgress(0.0);

        }

        public event EventHandler<ProgressEventArgs> UpdateProgress;

        protected virtual void OnUpdateProgress(double progress)
        {
            var handler = UpdateProgress;
            if (handler != null)
            {
                handler(this, new ProgressEventArgs(progress));
            }
        }

        public static void printStats(List<GlycoPSM> psms, string pathOut, List<double> glycanMasses)
        {
            // Make glyco peptide mass bins
            List<double> glycoPeptideMasses = new List<double>();
            //List<double> glycanMasses = new List<double>();
            List<double> peptideBackboneMasses = new List<double>();

            double maxGlycoPeptideMass = 0;
            double maxGlycanMass = 0;
            double maxPeptideBackboneMass = 0;

            foreach(GlycoPSM psm in psms)
            {
                glycoPeptideMasses.Add(psm.calcMH);
                //glycanMasses.Add(psm.peptide.GetModifications());
                peptideBackboneMasses.Add(psm.peptideMonoMass);

                if (psm.calcMH > maxGlycoPeptideMass)
                {
                    maxGlycoPeptideMass = psm.calcMH;
                }

                if(psm.peptideMonoMass > maxPeptideBackboneMass)
                {
                    maxPeptideBackboneMass = psm.peptideMonoMass;
                }
            }

            List<string> glycoPeptideMassBins = getBinCounts(glycoPeptideMasses, maxGlycoPeptideMass, 50);
            List<string> peptideBackboneMassBins = getBinCounts(peptideBackboneMasses, maxPeptideBackboneMass, 50);
            List<string> glycanMassBins = getBinCounts(glycanMasses, glycanMasses.Max(x => x), 200);


            // Make glycan mass bins
            

            StreamWriter writer = new StreamWriter(pathOut);
            writer.WriteLine("GlycoPeptide Mass,Count,,Peptide Backbone Mass,Count,,Glycan Mass,Count");


            bool stillOpen = true;
            int i = 0;

            while (stillOpen)
            {
                string nextLine = "";

                if(i < glycoPeptideMassBins.Count)
                {
                    nextLine += glycoPeptideMassBins[i] + ",,";
                }
                else
                {
                    nextLine += ",,,";
                }

                if (i < peptideBackboneMassBins.Count)
                {
                    nextLine += peptideBackboneMassBins[i] + ",,";
                }
                else
                {
                    nextLine += ",,,";
                }

                if(i < glycanMassBins.Count)
                {
                    nextLine += glycanMassBins[i] + ",,";
                }
                else
                {
                    nextLine += ",,,";
                }

                writer.WriteLine(nextLine);

                i++;

                if(i >= glycoPeptideMassBins.Count && i >= peptideBackboneMassBins.Count && i >= glycanMassBins.Count)
                {
                    stillOpen = false;
                }
            }

            writer.Close();

        }

        public static List<string> getBinCounts (List<double> numbers, double max, int binWidth)
        {
            //double min = 0;
            List<string> bins = new List<string>();

            for (int i = 0; i < max / binWidth + 5; i++)
            {
                int binCount = 0;

                foreach(double number in numbers)
                {
                    if(number > (i * binWidth) && number < (binWidth * (i + 1)))
                    {
                        binCount++;
                    }
                }
                
                bins.Add(i * binWidth + (binWidth / 2) + "," + binCount);

            }
           
            return bins;
        }

        public static string GetModName(string mod)
        {
            string[] modParse1 = mod.Split('(');
            string[] modParse2 = modParse1[1].Split(')');
            string[] modParse3 = modParse2[0].Split('/');
            string[] modParse4 = modParse3[0].Split(' ');

            string modName = modParse4[0];

            return modName;
        }

        public static double GetModMass(string mod)
        {
            string[] modParse1 = mod.Split('(');
            string[] modParse2 = modParse1[1].Split(')');
            string[] modParse3 = modParse2[0].Split('/');

            string massSubstring = modParse3[1].Substring(1);            

            double modMass = Convert.ToDouble(massSubstring);

            return modMass;
        }

        public static int GetModPosition(string mod)
        {
            string[] modParse1 = mod.Split('(');

            string modPositionString = modParse1[0].Substring(1);

            int modPosition = Convert.ToInt32(modPositionString);

            return modPosition;
        }

        public static double CalculateFDR(List<GlycoPSM> psmList)
        {
            double targetCount = 0;
            double decoyCount = 0;
            foreach (GlycoPSM psm in psmList)
            {
                if (psm.proteinName.Contains(">Reverse") || psm.proteinName.Contains("DECOY"))
                    decoyCount++;
                else
                    targetCount++;
            }

            double fdr = (decoyCount / targetCount) * 100.000000;

            return fdr;
        }

        public static double GetDeltaOfClosestRT(List<double> list)
        {
            double delta = 0;

            if (list.Count > 1)
            {
                List<double> sortedList = list.OrderBy(x => x).ToList();

                delta = sortedList[1] - sortedList[0];
                for (int i = 1; i < sortedList.Count - 1; i++)
                {
                    delta = Math.Min(delta, sortedList[i] - sortedList[i - 1]);
                }


            }

            return delta;
        }

        public static double GetDeltaOfFurthestRT(List<double> list)
        {
            double delta = 0;

            if (list.Count > 1)
            {
                List<double> sortedList = list.OrderBy(x => x).ToList();

                delta = sortedList[1] - sortedList[0];
                for (int i = 1; i < sortedList.Count; i++)
                {
                    delta = Math.Max(delta, sortedList[i] - sortedList[i - 1]);
                }


            }

            return delta;
        }

        public static double GetAvgDeltaRT(List<double> list)
        {
            double averageDeltaRT = 0;
            List<double> deltaList = new List<double>();

            if (list.Count > 1)
            {
                List<double> sortedList = list.OrderBy(x => x).ToList();

                double delta = sortedList[1] - sortedList[0];
                for (int i = 1; i < sortedList.Count; i++)
                {
                    delta = sortedList[i] - sortedList[i - 1];
                    deltaList.Add(delta);
                }
            }

            averageDeltaRT = Statistics.Mean(deltaList);

            return averageDeltaRT;
        }

        public static double GetStdDevDeltaRT(List<double> list)
        {
            double stdDevDeltaRT = 0;
            List<double> deltaList = new List<double>();

            if (list.Count > 1)
            {
                List<double> sortedList = list.OrderBy(x => x).ToList();

                double delta = sortedList[1] - sortedList[0];
                for (int i = 1; i < sortedList.Count; i++)
                {
                    delta = sortedList[i] - sortedList[i - 1];
                    deltaList.Add(delta);
                }
            }

            stdDevDeltaRT = Statistics.StandardDeviation(deltaList);

            return stdDevDeltaRT;
        }   

    }
}
