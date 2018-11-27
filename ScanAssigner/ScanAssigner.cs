using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using CSMSL;
using CSMSL.IO.Thermo;
using CSMSL.Proteomics;
using System.IO;
using LumenWorks.Framework.IO.Csv;

namespace ScanAssigner
{
    class ScanAssigner
    {
        public List<string> rawFilePaths;
        public List<string> byonicResultsPaths;
        public string outputFolderPath;
        public string uniprotGlycoDBFile;
        

        public ScanAssigner(List<string> rawFilePaths, List<string> byonicResultsPaths, string outputFolderPath, string uniprotGlycoDBFile)
        {
            this.rawFilePaths = rawFilePaths;
            this.byonicResultsPaths = byonicResultsPaths;
            this.outputFolderPath = outputFolderPath;
            this.uniprotGlycoDBFile = uniprotGlycoDBFile;
        }

        public void Start()
        {
            foreach (string file in byonicResultsPaths)
            {
                string rawFileToProcess = "";
                string resultsFile = Path.GetFileName(file);
                string namePart = resultsFile.Split('.')[0];

                foreach(string rawFile in rawFilePaths)
                {
                    if (rawFile.Contains(namePart + ".raw"))
                    {
                        rawFileToProcess = rawFile;
                    }                        
                }

                if (!rawFileToProcess.Equals(""))
                {
                    onHighlightListItems(file, rawFileToProcess);
                    crunch(rawFileToProcess, file, outputFolderPath);
                }
                else
                {
                    //MessageBox.Show("Cannot find raw file for " + file);
                }
            }
        }

        public void crunch(string rawFilePath, string byonicResultsPath, string outputFolderPath)
        {
            List<PSM> allPsmsQuant = new List<PSM>();
            OnUpdateProgress(0.25);

            ThermoRawFile rawFile = new ThermoRawFile(@rawFilePath);
            rawFile.Open();

            String file = Path.GetFileNameWithoutExtension(byonicResultsPath);

            StreamReader byonicFile = new StreamReader(byonicResultsPath);
            StreamWriter outputSummary = new StreamWriter(@outputFolderPath + "\\" + file + "_summary.csv");
            StreamWriter outputList = new StreamWriter(@outputFolderPath + "\\" + file + "_list.txt");


            //Write Headers
            outputList.WriteLine("Prot.Rank\tSequence\tPeptideParseFriendly\tPeptide\tPos.\tMods(variable)\tGlycans\t"+
                                 "PEP2D\tPEP1D\t|Log Prob|\tScore\tDeltaScore\tDelta Mod.Score\tCharge\tObs.m/z\t"+
                                 "Precursor Theoretical m/z (Th)\tppmerr.\tObs.MH\tCalc.MH\tCleavage\tGlycansPos.\tProteinName" +
                                 "\tProt.Id\tScanTime\tSpectrum Number\tMods(fixed)\tFDR2D\tFDR1D\t"+
                                 "FDR uniq.2D\tFDR uniq.1D\tq-value2D\tq-value1D\tisGlycoPeptide\tmodsPassedCheck" +
                                 "\tpositionPassedCheck\tfragmentation\tMasterScan");

            List<int> HCDscans = new List<int>();
            List<int> ETDscans = new List<int>();
            int totalMSMSscans = 0;

            //these will have decoys and targets, so only use is for FDR calc (FDR calc will be at PSM level)
            List<PSM> allPSMs = new List<PSM>();
            List<PSM> hcdPSMs = new List<PSM>();
            List<PSM> etdPSMs = new List<PSM>();
            List<PSM> allGlycoPSMs = new List<PSM>();
            List<PSM> hcdGlycoPSMs = new List<PSM>();
            List<PSM> etdGlycoPSMs = new List<PSM>();

            //use this to print file to use in next program
            //List<PSM> masterListPSMs = new List<PSM>();

            //these will have no decoys, so use for counts and other processing
            Dictionary<string, List<PSM>> allPSMsNODECOY = new Dictionary<string, List<PSM>>();
            Dictionary<string, List<PSM>> hcdPSMsNODECOY = new Dictionary<string, List<PSM>>();
            Dictionary<string, List<PSM>> etdPSMsNODECOY = new Dictionary<string, List<PSM>>();
            Dictionary<string, List<PSM>> allGlycoPSMsNODECOY = new Dictionary<string, List<PSM>>();
            Dictionary<string, List<PSM>> hcdGlycoPSMsNODECOY = new Dictionary<string, List<PSM>>();
            Dictionary<string, List<PSM>> etdGlycoPSMsNODECOY = new Dictionary<string, List<PSM>>();


            for (int i = rawFile.FirstSpectrumNumber; i <= rawFile.LastSpectrumNumber; i++)
            {
                if (rawFile.GetMsnOrder(i) == 2)
                {
                    if (rawFile.GetDissociationType(i).ToString().Equals("HCD"))
                    {
                        HCDscans.Add(i);
                    }
                    if (rawFile.GetDissociationType(i).ToString().Equals("ETD"))
                    {
                        ETDscans.Add(i);
                    }
                    totalMSMSscans++;
                }
            }

            

            using (var csv = new CsvReader(byonicFile, true, '\t'))
            {
                while (csv.ReadNextRecord())
                {
                    try
                    {
                        List<int> glycanPositionsList = new List<int>();

                        string PID = "NA";//csv["PID"];
                        string protRank = csv["ProtRank"];
                        string byonicSequence = csv["Sequence"];
                        string peptidesToBeParsed = csv["PeptideParseFriendly"];

                        // NEW
                        string aminoAcidSequence = peptidesToBeParsed.Split(',')[0];
                        string byonicIntensity = "NA"; //csv["Intensity"];

                        int peptideStartPosition = int.Parse(csv["Pos."]);
                        string modsToBeParsed = csv["Mods(variable)"];
                        string glycansToBeParsed = csv["Glycans"];
                        double PEP2D = double.Parse(csv["PEP2D"]);
                        double PEP1D = double.Parse(csv["PEP1D"]);
                        //double logProb = double.Parse(csv["|Log Prob|"]);
                        double logProb = -1;
                        double score = double.Parse(csv["Score"]);
                        double deltaScore = double.Parse(csv["DeltaScore"]);
                        double deltaModScore = double.Parse(csv["Delta Mod.Score"]);
                        int charge = int.Parse(csv["z"]);
                        double mzObs = double.Parse(csv["Obs.m/z"]);
                        double mzCalc = double.Parse(csv["Calc.m/z"]);
                        //double ppmError = double.Parse(csv["ppmerr."]);
                        string ppmError = csv["ppmerr."];
                        double obsMH = double.Parse(csv["Obs.MH"]);
                        double calcMH = double.Parse(csv["Calc.MH"]);
                        string cleavage = csv["Cleavage"];
                        string glycanPositions = csv["GlycansPos."];
                        string proteinName = csv["ProteinName"];
                        int protID = int.Parse(csv["Prot.Id"]);
                        double scanTime = double.Parse(csv["ScanTime"]);                        
                        int scanNumber = int.Parse(csv["Scan #"]);
                        string modsFixed = csv["Mods(fixed)"];
                        double FDR2D = double.Parse(csv["FDR2D"]);
                        double FDR1D = double.Parse(csv["FDR1D"]);
                        double FDR2Dunique = double.Parse(csv["FDR uniq.2D"]);
                        double FDR1Dunique = double.Parse(csv["FDR uniq.1D"]);
                        double qvalue2D = double.Parse(csv["q-value2D"]);
                        double qvalue1D = double.Parse(csv["q-value1D"]);

                        // Flags to be used
                        bool isGlycoPeptide = false;
                        bool modsPassedCheck = false;
                        bool positionsPassedCheck = false;

                        if (modsToBeParsed.Contains("Glycan"))
                        {
                            isGlycoPeptide = true;
                        }

                        string fragmentation = rawFile.GetDissociationType(scanNumber).ToString();

                        //ThermoSpectrum spectrum = rawFile.GetSpectrum(scanNumber);

                        string[] parsedPeptide = peptidesToBeParsed.Split(',');
                        Peptide peptide = new Peptide(parsedPeptide[0]);

                        //string[] parsedSequenceWithMods = sequence.Split('.');
                        //string uniqueSequence = parsedSequenceWithMods[1];

                        string uniqueSequence = byonicSequence.Substring(2, byonicSequence.Length - 4);
                        //Console.WriteLine(uniqueSequence);
                        //Console.WriteLine(newSequence);
                        //Console.ReadKey();
                        //add to master list here

                        if (!String.IsNullOrEmpty(modsToBeParsed))
                        {
                            string[] modsParsedArrary = modsToBeParsed.Split(';');
                            for (int i = 0; i < modsParsedArrary.Length; i++)
                            {
                                string mod = modsParsedArrary[i];
                                //Console.WriteLine(mod);
                                //Console.ReadKey();
                                if (mod[0].Equals(' '))
                                {
                                    mod = mod.Substring(1);
                                }

                                if (mod.Contains("NGlycan") && !mod[0].Equals('N'))
                                {
                                    modsPassedCheck = false;
                                }

                                if (mod.Contains("NGlycan") && mod[0].Equals('N'))
                                {
                                    modsPassedCheck = true;
                                }

                                if (mod.Contains("OGlycan") && mod[0].Equals('S'))
                                {
                                    modsPassedCheck = true;
                                }
                                if (mod.Contains("OGlycan") && mod[0].Equals('T'))
                                {
                                    modsPassedCheck = true;
                                }
                            }
                        }

                        if (!String.IsNullOrEmpty(glycanPositions))
                        {
                            string[] glycansPosParsedArray = glycanPositions.Split(';');
                            for (int i = 0; i < glycansPosParsedArray.Length; i++)
                            {
                                int glycanPos = Convert.ToInt32(glycansPosParsedArray[i]);
                                glycanPositionsList.Add(glycanPos);
                            }
                        }

                        OnUpdateProgress(0.5);
                        if (glycanPositionsList.Count != glycanPositionsList.Distinct().Count())
                        {
                            positionsPassedCheck = false;
                        }
                        else
                        {
                            positionsPassedCheck = true;
                        }

                        int masterScan = rawFile.GetParentSpectrumNumber(scanNumber);

                        PSM psm = new PSM(peptide, PID, protRank, peptidesToBeParsed, peptideStartPosition, modsToBeParsed, glycansToBeParsed, 
                                          PEP2D, PEP1D, logProb, score, deltaScore, deltaModScore, charge, mzObs, mzCalc, ppmError,
                                          obsMH, calcMH, cleavage, glycanPositions, proteinName, protID, scanTime, scanNumber, modsFixed, FDR2D, FDR1D, 
                                          FDR2Dunique, FDR1Dunique, qvalue2D, qvalue1D, isGlycoPeptide, fragmentation, byonicSequence, 
                                          modsPassedCheck, masterScan, positionsPassedCheck, byonicIntensity);//, spectrum);


                        if (peptide.Length > 4 && FDR2D <= 0.01)
                        {
                            allPSMs.Add(psm);
                            
                            if (fragmentation.Equals("HCD"))
                            {
                                hcdPSMs.Add(psm);
                            }
                            if (fragmentation.Equals("ETD"))
                            {
                                etdPSMs.Add(psm);
                            }

                            if (isGlycoPeptide)
                            {
                                allGlycoPSMs.Add(psm);
                                if (fragmentation.Equals("HCD"))
                                {
                                    hcdGlycoPSMs.Add(psm);
                                }
                                if (fragmentation.Equals("ETD"))
                                {
                                    etdGlycoPSMs.Add(psm);
                                }
                            }

                            if (!proteinName.Contains("DECOY") && !proteinName.Contains("Reverse") && modsPassedCheck && positionsPassedCheck)
                            {
                                if (!allPSMsNODECOY.ContainsKey(uniqueSequence))
                                {
                                    allPSMsNODECOY.Add(uniqueSequence, new List<PSM>());
                                }
                                allPSMsNODECOY[uniqueSequence].Add(psm);

                                if (fragmentation.Equals("HCD"))
                                {
                                    if (!hcdPSMsNODECOY.ContainsKey(uniqueSequence))
                                    {
                                        hcdPSMsNODECOY.Add(uniqueSequence, new List<PSM>());
                                    }
                                    hcdPSMsNODECOY[uniqueSequence].Add(psm);
                                }

                                if (fragmentation.Equals("ETD"))
                                {
                                    if (!etdPSMsNODECOY.ContainsKey(uniqueSequence))
                                    {
                                        etdPSMsNODECOY.Add(uniqueSequence, new List<PSM>());
                                    }
                                    etdPSMsNODECOY[uniqueSequence].Add(psm);
                                }

                                if (isGlycoPeptide)
                                {
                                    if (!allGlycoPSMsNODECOY.ContainsKey(uniqueSequence))
                                    {
                                        allGlycoPSMsNODECOY.Add(uniqueSequence, new List<PSM>());
                                    }
                                    allGlycoPSMsNODECOY[uniqueSequence].Add(psm);

                                    if (fragmentation.Equals("HCD"))
                                    {
                                        if (!hcdGlycoPSMsNODECOY.ContainsKey(uniqueSequence))
                                        {
                                            hcdGlycoPSMsNODECOY.Add(uniqueSequence, new List<PSM>());
                                        }
                                        hcdGlycoPSMsNODECOY[uniqueSequence].Add(psm);
                                    }

                                    if (fragmentation.Equals("ETD"))
                                    {
                                        if (!etdGlycoPSMsNODECOY.ContainsKey(uniqueSequence))
                                        {
                                            etdGlycoPSMsNODECOY.Add(uniqueSequence, new List<PSM>());
                                        }
                                        etdGlycoPSMsNODECOY[uniqueSequence].Add(psm);
                                    }
                                }
                            }
                            
                        }

                        /**
                        psm.Peak126 = GetPeak(psm, 126.055, rawFile);
                        psm.Peak138 = GetPeak(psm, 138.055, rawFile);
                        psm.Peak144 = GetPeak(psm, 144.065, rawFile);
                        psm.Peak168 = GetPeak(psm, 168.066, rawFile);
                        psm.Peak186 = GetPeak(psm, 186.076, rawFile);
                        psm.Peak204 = GetPeak(psm, 204.087, rawFile);
                        psm.Peak274 = GetPeak(psm, 274.092, rawFile);
                        psm.Peak292 = GetPeak(psm, 292.103, rawFile);
                        psm.Peak366 = GetPeak(psm, 366.140, rawFile);
                        **/

                        double GlcNAcToGalNAcratio = 0;

                        /*
                        outputList.Write("PID\tProt.Rank\tSequence\tPeptideParseFriendly\tPos.\tMods(variable)\tGlycans\tPEP2D\tPEP1D\t|Log Prob|\tScore\tDeltaScore\tDelta Mod.Score\tz\tObs.m/z\t");
                        outputList.Write("Calc.m/z\tppmerr.\tObs.m/z\tCalc.m/z\tCleavage\tGlycansPos.\tProteinName\tProt.Id\tScanTime\tScan #\tMods(fixed)\tFDR2D\tFDR1D\t");
                        outputList.WriteLine("FDR uniq.2D\tFDR uniq.1D\tq-value2D\tq-value1D\tisGlycoPeptide\tmodsPassedCheck\tpositionPassedCheck\tfragmentation\tPeak126\tPeak138\tPeak144\tPeak168\tPeak186\tPeak204\tPeak274\tPeak292\tPeak366");
                        */

                        outputList.WriteLine(protRank + "\t" + byonicSequence + "\t" + peptidesToBeParsed + "\t" + aminoAcidSequence + "\t" +
                                             peptideStartPosition + "\t" + modsToBeParsed + "\t" + glycansToBeParsed + "\t" + PEP2D + "\t" + PEP1D + "\t" +
                                             logProb + "\t" + score + "\t" + deltaScore + "\t" + deltaModScore + "\t" + charge + "\t" + mzObs + "\t" + mzCalc +
                                             "\t" + ppmError + "\t" + obsMH + "\t" + calcMH + "\t" + cleavage + "\t" + glycanPositions + "\t" + proteinName + "\t" +
                                             protID + "\t" + scanTime + "\t" + scanNumber + "\t" + modsFixed + "\t" + FDR2D + "\t" + FDR1D + "\t" + FDR2Dunique +
                                             "\t" + FDR1Dunique + "\t" + qvalue2D + "\t" + qvalue1D + "\t" + isGlycoPeptide + "\t" + modsPassedCheck + "\t" +
                                             positionsPassedCheck + "\t" + fragmentation + "\t" + masterScan);

                    } catch(Exception e)
                    {
                        Console.WriteLine();
                    }                    
                }
            }

            //24 7,Console.ReadKey();
            //For Nick's LFQ
            try
            {
                LFQProcessor lfqs = new LFQProcessor(@rawFilePath);
                allPsmsQuant = lfqs.crunch(allPSMs);
            }
            catch (Exception e)
            {
                System.Windows.Forms.MessageBox.Show(e.ToString());
            }

            // Start LFQ Writer
            StreamWriter quantWriter = new StreamWriter(@outputFolderPath + "\\" + file + "_list_Quant.txt");

            quantWriter.WriteLine("Prot.Rank\tSequence\tPeptideParseFriendly\tPeptide\tPos.\tMods(variable)\tGlycans"+
                                  "\tPEP2D\tPEP1D\t|Log Prob|\tScore\tDeltaScore\tDelta Mod.Score\tCharge\tObs.m/z\t" +
                                  "Precursor Theoretical m/z (Th)\tppmerr.\tObs.MH\tCalc.MH\tCleavage\tGlycansPos.\tProteinName"+
                                  "\tProt.Id\tScanTime\tSpectrum Number\tMods(fixed)\tFDR2D\tFDR1D\t" +
                                  "FDR uniq.2D\tFDR uniq.1D\tq-value2D\tq-value1D\tisGlycoPeptide\tmodsPassedCheck"+
                                  "\tpositionPassedCheck\tfragmentation\tMasterScan\tLFQ Intensity");

            foreach (PSM psmQ in allPsmsQuant)
            {
                quantWriter.WriteLine(psmQ.protRank + "\t" + psmQ.byonicSequence + "\t" + psmQ.peptidesToBeParsed + "\t" + psmQ.peptide.Sequence + "\t" +
                    psmQ.peptideStartPosition + "\t" + psmQ.modsToBeParsed + "\t" + psmQ.glycansToBeParsed + "\t" + psmQ.PEP2D + "\t" + psmQ.PEP1D + "\t" + psmQ.logProb +
                    "\t" + psmQ.score + "\t" + psmQ.deltaScore + "\t" + psmQ.deltaModScore + "\t" + psmQ.charge + "\t" + psmQ.mzObs + "\t" + psmQ.mzCalc + "\t" + psmQ.ppmError +
                    "\t" + psmQ.obsMH + "\t" + psmQ.calcMH + "\t" + psmQ.cleavage + "\t" + psmQ.glycanPositions + "\t" + psmQ.proteinName + "\t" + psmQ.protID + "\t" + psmQ.scanTime +
                    "\t" + psmQ.scanNumber + "\t" + psmQ.modsFixed + "\t" + psmQ.FDR2D + "\t" + psmQ.FDR1D + "\t" + psmQ.FDR2Dunique + "\t" + psmQ.FDR1Dunique + "\t" +
                    psmQ.qvalue2D + "\t" + psmQ.qvalue1D + "\t" + psmQ.isGlycoPeptide + "\t" + psmQ.modsPassedCheck + "\t" + psmQ.positionPassedCheck + "\t" + psmQ.fragmentation
                    + "\t" + psmQ.masterScan + "\t" + psmQ.intensity);
            }

            quantWriter.Close();
            // End LFQ Writer


            OnUpdateProgress(0.75);
            double allPSMsFDR = CalculateFDR(allPSMs);
            double hcdPSMsFDR = CalculateFDR(hcdPSMs);
            double etdPSMsFDR = CalculateFDR(etdPSMs);
            double allGlycoPSMsFDR = CalculateFDR(allGlycoPSMs);
            double hcdGlycoPSMsFDR = CalculateFDR(hcdGlycoPSMs);
            double etdGlycoPSMsFDR = CalculateFDR(etdGlycoPSMs);

            int allPSMCount = GetPSMcount(allPSMsNODECOY);
            int hcdPSMCount = GetPSMcount(hcdPSMsNODECOY);
            int etdPSMCount = GetPSMcount(etdPSMsNODECOY);
            int allGlycoPSMCount = GetPSMcount(allGlycoPSMsNODECOY);
            int hcdGlycoPSMCount = GetPSMcount(hcdGlycoPSMsNODECOY);
            int etdGlycoPSMCount = GetPSMcount(etdGlycoPSMsNODECOY);

            int allPeptideCount = allPSMsNODECOY.Count;
            int hcdPeptideCount = hcdPSMsNODECOY.Count;
            int etdPeptideCount = etdPSMsNODECOY.Count;
            int allGlycoPeptideCount = allGlycoPSMsNODECOY.Count;
            int hcdGlycoPeptideCount = hcdGlycoPSMsNODECOY.Count;
            int etdGlycoPeptideCount = etdGlycoPSMsNODECOY.Count;

            int hcdWithHexNAcOxonium = GetSpectraCountOxonium(HCDscans, 204.087, 0, rawFile);
            int hcdWithHexNAcHexOxonium = GetSpectraCountOxonium(HCDscans, 366.14, 0, rawFile);
            int hcdWithHexNAcHex2Oxonium = GetSpectraCountOxonium(HCDscans, 528.19, 0, rawFile);
            int hcdWithHexNAcHex3Oxonium = GetSpectraCountOxonium(HCDscans, 690.245, 0, rawFile);
            int hcdWithHexNAc2Hex3Oxonium = GetSpectraCountOxonium(HCDscans, 893.325, 0, rawFile);
            int hcdWithHexNAcHexNeuAcOxonium = GetSpectraCountOxonium(HCDscans, 657.235, 0, rawFile);
            int hcdWithHexoseOxonium = GetSpectraCountOxonium(HCDscans, 163.060, 0, rawFile);
            int hcdWithFucOxonium = GetSpectraCountOxonium(HCDscans, 147.065, 0, rawFile);
            int hcdWithNeuAcOxonium = GetSpectraCountOxonium(HCDscans, 274.092, 292.103, rawFile);
            int hcdWithNeuAcAcetylOxonium = GetSpectraCountOxonium(HCDscans, 316.103, 334.113, rawFile);
            int hcdWithHexNAcNeuAcOxonium = GetSpectraCountOxonium(HCDscans, 495.182, 0, rawFile);
            int hcdWithNeuGcOxonium = GetSpectraCountOxonium(HCDscans, 290.087, 308.098, rawFile);
            int hcdWithNeuGcAcetylOxonium = GetSpectraCountOxonium(HCDscans, 332.098, 350.108, rawFile);
            int hcdWithPhosphoHexOxonium = GetSpectraCountOxonium(HCDscans, 243.0264, 0, rawFile);
            int etdWithHexNAcOxonium = GetSpectraCountOxonium(ETDscans, 204.087, 0, rawFile);
            int etdWithHexNAcHexOxonium = GetSpectraCountOxonium(ETDscans, 366.14, 0, rawFile);
            int etdWithHexNAcHex2Oxonium = GetSpectraCountOxonium(ETDscans, 528.19, 0, rawFile);
            int etdWithHexNAcHex3Oxonium = GetSpectraCountOxonium(ETDscans, 690.245, 0, rawFile);
            int etdWithHexNAc2Hex3Oxonium = GetSpectraCountOxonium(ETDscans, 893.325, 0, rawFile);
            int etdWithHexNAcHexNeuAcOxonium = GetSpectraCountOxonium(ETDscans, 657.235, 0, rawFile);
            int etdWithHexoseOxonium = GetSpectraCountOxonium(ETDscans, 163.060, 0, rawFile);
            int etdWithFucOxonium = GetSpectraCountOxonium(ETDscans, 147.065, 0, rawFile);
            int etdWithNeuAcOxonium = GetSpectraCountOxonium(ETDscans, 274.092, 292.103, rawFile);
            int etdWithNeuAcAcetylOxonium = GetSpectraCountOxonium(ETDscans, 316.103, 334.113, rawFile);
            int etdWithHexNAcNeuAcOxonium = GetSpectraCountOxonium(ETDscans, 495.182, 0, rawFile);
            int etdWithNeuGcOxonium = GetSpectraCountOxonium(ETDscans, 290.087, 308.098, rawFile);
            int etdWithNeuGcAcetylOxonium = GetSpectraCountOxonium(ETDscans, 332.098, 350.108, rawFile);
            int etdWithPhosphoHexOxonium = GetSpectraCountOxonium(ETDscans, 243.0264, 0, rawFile);
            OnUpdateProgress(1.0);

            outputSummary.WriteLine(",Total,HCD,ETD");
            outputSummary.WriteLine("MS/MS Scans," + totalMSMSscans + "," + HCDscans.Count + "," + ETDscans.Count);
            outputSummary.WriteLine("PSMs (all)," + allPSMCount + "," + hcdPSMCount + "," + etdPSMCount);
            outputSummary.WriteLine("Unique Peptides (all)," + allPeptideCount + "," + hcdPeptideCount + "," + etdPeptideCount);
            outputSummary.WriteLine("FDR PSM level," + allPSMsFDR + "," + hcdPSMsFDR + "," + etdPSMsFDR);
            outputSummary.WriteLine("GlycoPSMs," + allGlycoPSMCount + "," + hcdGlycoPSMCount + "," + etdGlycoPSMCount);
            outputSummary.WriteLine("Unique GlycoPeptides," + allGlycoPeptideCount + "," + hcdGlycoPeptideCount + "," + etdGlycoPeptideCount);
            outputSummary.WriteLine("FDR GlycoPSM level," + allGlycoPSMsFDR + "," + hcdGlycoPSMsFDR + "," + etdGlycoPSMsFDR);
            outputSummary.WriteLine("Oxonium Ions /////,/////////,/////////,/////////");
            outputSummary.WriteLine("HexNAc (204.09),," + hcdWithHexNAcOxonium + "," + etdWithHexNAcOxonium);
            outputSummary.WriteLine("HexNAcHex (366.14),," + hcdWithHexNAcHexOxonium + "," + etdWithHexNAcHexOxonium);
            outputSummary.WriteLine("HexNAcHex(2) (528.19),," + hcdWithHexNAcHex2Oxonium + "," + etdWithHexNAcHex2Oxonium);
            outputSummary.WriteLine("HexNAcHex(3) (690.245),," + hcdWithHexNAcHex3Oxonium + "," + etdWithHexNAcHex3Oxonium);
            outputSummary.WriteLine("HexNAc(2)Hex(3) (893.325),," + hcdWithHexNAc2Hex3Oxonium + "," + etdWithHexNAc2Hex3Oxonium);
            outputSummary.WriteLine("HexNAcHexNeuAc (657.23),," + hcdWithHexNAcHexNeuAcOxonium + "," + etdWithHexNAcHexNeuAcOxonium);
            outputSummary.WriteLine("Hexose (163.06),," + hcdWithHexoseOxonium + "," + etdWithHexoseOxonium);
            outputSummary.WriteLine("Fucose (147.065),," + hcdWithFucOxonium + "," + etdWithFucOxonium);
            outputSummary.WriteLine("NeuAc (272.09/292.103),," + hcdWithNeuAcOxonium + "," + etdWithNeuAcOxonium);
            outputSummary.WriteLine("NeuAcAcetyl (316.103/334.113),," + hcdWithNeuAcAcetylOxonium + "," + etdWithNeuAcAcetylOxonium);
            outputSummary.WriteLine("HexNAcNeuAc (495.1382),," + hcdWithHexNAcNeuAcOxonium + "," + etdWithHexNAcNeuAcOxonium);
            outputSummary.WriteLine("NeuGc (290.087/308.098),," + hcdWithNeuGcOxonium + "," + etdWithNeuGcOxonium);
            outputSummary.WriteLine("NeuGcAcetyl (332.098/350.108),," + hcdWithNeuGcAcetylOxonium + "," + etdWithNeuGcAcetylOxonium);
            outputSummary.WriteLine("PhosphoHex (243.0264),," + hcdWithPhosphoHexOxonium + "," + etdWithPhosphoHexOxonium);

            outputSummary.Close();
            outputList.Close();

            var glycoFragmentFinder = new glycoFragmentFinder(rawFile, "1", @outputFolderPath + "\\" + file + "_list.txt", uniprotGlycoDBFile, outputFolderPath);
            glycoFragmentFinder.crunch();


            rawFile.Dispose();

            OnUpdateProgress(0.0);
        }

        private static double CalculateFDR(List<PSM> psmList)
        {
            double targetCount = 0;
            double decoyCount = 0;
            foreach (PSM psm in psmList)
            {
                if (psm.proteinName.Contains(">Reverse"))
                    decoyCount++;
                else
                    targetCount++;
            }

            double fdr = (decoyCount / targetCount) * 100.000000;

            //List<PSM> psmListFDRorder = psmList.OrderBy(x => x.FDR2D).ToList();

            return fdr;
        }

        private static int GetPSMcount(Dictionary<string, List<PSM>> dictionary)
        {
            int count = 0;
            foreach (KeyValuePair<string, List<PSM>> entry in dictionary)
            {
                foreach (PSM psm in entry.Value)
                {
                    count++;
                }
            }
            return count;
        }

        public static int GetSpectraCountOxonium(List<int> spectrumList, double oxoniumIon1, double oxoniumIon2, ThermoRawFile rawFile)
        {
            int count = 0;


            foreach (int scan in spectrumList)
            {

                ThermoSpectrum spectrum = rawFile.GetSpectrum(scan);

                bool oxoniumIsPresent = false;
                DoubleRange rangeOxonium1 = DoubleRange.FromPPM(oxoniumIon1, 20);
                List<ThermoMzPeak> peaks;
                if (spectrum.TryGetPeaks(rangeOxonium1, out peaks))
                {
                    peaks = peaks.OrderBy(x => x.SignalToNoise).ToList();
                }

                if (peaks.Count > 0)
                {
                    if (peaks[0].SignalToNoise > 3)
                    {
                        oxoniumIsPresent = true;
                    }
                }

                if (oxoniumIon2 > 0)
                {
                    DoubleRange rangeOxonium2 = DoubleRange.FromPPM(oxoniumIon2, 20);
                    List<ThermoMzPeak> peaks2;
                    if (spectrum.TryGetPeaks(rangeOxonium2, out peaks2))
                    {
                        peaks2 = peaks2.OrderBy(x => x.SignalToNoise).ToList();
                    }

                    if (peaks2.Count > 0)
                    {
                        if (peaks2[0].SignalToNoise > 3)
                        {
                            oxoniumIsPresent = true;
                        }
                    }
                }

                if (oxoniumIsPresent)
                {
                    count++;
                }
            }


            return count;
        }

        private static ThermoMzPeak GetPeak(PSM psm, double oxoniumIon, ThermoRawFile rawFile)
        {
            DoubleRange rangeOxonium = DoubleRange.FromPPM(oxoniumIon, 20);
            List<ThermoMzPeak> peaks;

            ThermoSpectrum spec = rawFile.GetSpectrum(psm.scanNumber);

            if (spec.TryGetPeaks(rangeOxonium, out peaks))
            {
                peaks = peaks.OrderBy(x => x.SignalToNoise).ToList();
            }

            double diff = double.MaxValue;
            ThermoMzPeak returnPeak = null;
            foreach (ThermoMzPeak peak in peaks)
            {
                var currDiff = Math.Abs(peak.MZ - oxoniumIon);
                if (currDiff < diff)
                {
                    diff = currDiff;
                    returnPeak = peak;
                }
            }
            return returnPeak;
        }

        public static double CheckPeak(ThermoMzPeak peak)
        {
            double intensity = 0;

            if (peak != null)
            {
                intensity = peak.Intensity;
            }

            return intensity;
        }

        private void scanAssignProgress_Click(object sender, EventArgs e)
        {

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

        public EventHandler<HighlightEventArgs> highlightListItems;

        protected virtual void onHighlightListItems( string resultFile, string rawFile)
        {
            var handler = highlightListItems;

            if(handler != null)
            {
                handler(this, new HighlightEventArgs(resultFile, rawFile));
            }
        }
    }
}
