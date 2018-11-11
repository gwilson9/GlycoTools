using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using CSMSL;
using CSMSL.IO.Thermo;
using CSMSL.Spectral;
using CSMSL.Util;
using LumenWorks.Framework.IO.Csv;
using CSMSL.Analysis;
using CSMSL.Proteomics;
using CSMSL.Chemistry;

namespace GlycoFragmentFinder
{
    class Program
    {
        static void Main(string[] args)
        {

            //Read in raw file
            ThermoRawFile rawFile = new ThermoRawFile(@"F:\Spring 2017\Glyco\2017_01_18_Data\Glyco EThcD\2017_01_19_MB1-4_Glyco_EThcD25_rep1_Frac12.raw");
            rawFile.Open();

            string fraction = "Rep1_Frac12";

            StreamReader inputFile = new StreamReader(@"F:\Spring 2017\Glyco\2017_01_18_Data\April 2017 Processing\ProtDB5 GlycoDB6 use\EThcD\EThcD_ByonicModMotifPresent_Combo_Frac12_list.txt");

            StreamReader UniprotGlycoDBfile = new StreamReader(@"F:\Spring 2017\Glyco\Algorithm Dev\Database Builder\Uniprot GlycoSites\MouseGlycoSites.csv");

            string outputPath = @"F:\Spring 2017\Glyco\2017_01_18_Data\April 2017 Processing\ProtDB5 GlycoDB6 use\EThcD\EThcD_ByonicModMotifPresent_Combo_Frac12_list";

            int scoreFilter = 50;

            StreamWriter outputSummary = new StreamWriter(outputPath + "_scoredSpectra_Filter" + scoreFilter + ".txt");
            StreamWriter outputEachFragment = new StreamWriter(outputPath + "_allFragments_Filter" + scoreFilter + ".txt");
            StreamWriter outputPeptideFragments = new StreamWriter(outputPath + "_peptideFragments_Filter" + scoreFilter + ".txt");
            StreamWriter outputGlycanFragments = new StreamWriter(outputPath + "_glycanFragments_Filter" + scoreFilter + ".txt");
            StreamWriter outputOxoniumFragments = new StreamWriter(outputPath + "_oxoniumFragments_Filter" + scoreFilter + ".txt");

            outputSummary.Write("Peptide\tSequenceOnly\tGlycan\tFraction\tScanNumber\tFragmentation\tMZ\tCharge\tPeptideMass\tGlycanMass\tScore\tDeltaModScore\t2DFDR\t2DPEP\tSeqCoverage\tSeqCovNoNeutLoss\tGlycanSeqCov\t");
            outputSummary.Write("NumberOfFragments\tNeutralLossFragments\tGlycanFragments\tOxoniumFragments\tFragWithGlycan\tNtermFragWithGlycan\tCtermFragWithGlycan\tTotalFragmentIntensity\tNoLossIntensity\t");
            outputSummary.WriteLine("GlycanFragIntensity\tOxoniumIntensity\tTIC\tFragmentList\tGlycanFragmentList\tOxoniumList");

            outputEachFragment.Write("Peptide\tSequenceOnly\tGlycan\tFraction\tScanNumber\tFragmentation\tMZ\tCharge\tPeptideMass\tGlycanMass\tScore\tDeltaModScore\t2DFDR\t2DPEP\tSeqCoverage\tSeqCovNoNeutLoss\tGlycanSeqCov\t");
            outputEachFragment.WriteLine("NumberOfFragments\tNeutralLossFragments\tGlycanFragments\tOxoniumFragments\tFragmentInfo");

            outputPeptideFragments.Write("Peptide\tSequenceOnly\tGlycan\tFraction\tScanNumber\tFragmentation\tMZ\tCharge\tPeptideMass\tGlycanMass\tScore\tDeltaModScore\t2DFDR\t2DPEP\tSeqCoverage\tSeqCovNoNeutLoss\tGlycanSeqCov\t");
            outputPeptideFragments.WriteLine("NumberOfFragments\tNeutralLossFragments\tFragmentInfo");

            outputGlycanFragments.Write("Peptide\tSequenceOnly\tGlycan\tFraction\tScanNumber\tFragmentation\tMZ\tCharge\tPeptideMass\tGlycanMass\tScore\tDeltaModScore\t2DFDR\t2DPEP\tSeqCoverage\tSeqCovNoNeutLoss\tGlycanSeqCov\t");
            outputGlycanFragments.WriteLine("GlycanFragments\tFragmentInfo");

            outputOxoniumFragments.Write("Peptide\tSequenceOnly\tGlycan\tFraction\tScanNumber\tFragmentation\tMZ\tCharge\tPeptideMass\tGlycanMass\tScore\tDeltaModScore\t2DFDR\t2DPEP\tSeqCoverage\tSeqCovNoNeutLoss\tGlycanSeqCov\t");
            outputOxoniumFragments.WriteLine("OxoniumFragments\tFragmentInfo");

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

            List<GlycoPSM> glycoPSMsLocalized = new List<GlycoPSM>();

            using (var csv = new CsvReader(inputFile, true, '\t'))
            {
                while (csv.ReadNextRecord())
                {

                    List<int> glycanPositionsList = new List<int>();
                    List<string> mods = new List<string>();
                    List<string> glycans = new List<string>();

                    string PID = csv["PID"];
                    string protRank = csv["Prot.Rank"];
                    string sequence = csv["Sequence"];
                    string peptidesToBeParsed = csv["PeptideParseFriendly"];
                    int peptideStartPosition = int.Parse(csv["Pos."]);
                    string modsToBeParsed = csv["Mods(variable)"];
                    string glycansToBeParsed = csv["Glycans"];
                    double PEP2D = double.Parse(csv["PEP2D"]);
                    double logProb = double.Parse(csv["|Log Prob|"]);
                    double score = double.Parse(csv["Score"]);
                    double deltaScore = double.Parse(csv["DeltaScore"]);
                    double deltaModScore = double.Parse(csv["Delta Mod.Score"]);
                    int charge = int.Parse(csv["z"]);
                    double mzObs = double.Parse(csv["Obs.m/z"]);
                    double ppmError = double.Parse(csv["ppmerr."]);
                    double obsMH = double.Parse(csv["Obs.MH"]);
                    string cleavage = csv["Cleavage"];
                    string glycanPositions = csv["GlycansPos."];
                    string proteinName = csv["ProteinName"];
                    int protID = int.Parse(csv["Prot.Id"]);
                    double scanTime = double.Parse(csv["ScanTime"]);
                    int scanNumber = int.Parse(csv["Scan #"]);
                    string modsFixed = csv["Mods(fixed)"];
                    double FDR2D = double.Parse(csv["FDR2D"]);
                    double FDR2Dunique = double.Parse(csv["FDR uniq.2D"]);
                    double qvalue2D = double.Parse(csv["q-value2D"]);
                    bool isGlycoPeptide = bool.Parse(csv["isGlycoPeptide"]);
                    string fragmentation = csv["fragmentation"];

                    bool seenWithHCD = false;
                    bool seenWithETD = false;
                    bool NXSmotif = false;
                    bool NXTmotif = false;
                    bool isLocalized = false;
                    bool Nlinked = false;
                    bool Olinked = false;

                    List<double> glycanMasses = new List<double>();

                    bool matchedToUniprot = false;
                    string uniprotEvidenceType = "None";
                    int uniprotEvidenceNumber = 0;

                    string[] parsedPeptide = peptidesToBeParsed.Split(',');
                    Peptide peptide = new Peptide(parsedPeptide[0]);
                    Peptide peptideNoGlycan = new Peptide(parsedPeptide[0]);                //added to enable ability to look for glycan (currenlty not used)
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



                    if (isGlycoPeptide && !proteinName.Contains("DECOY") && !proteinName.Contains("Reverse") && peptide.Length > 4 && FDR2D <= 0.01 && score >= scoreFilter)          //!!!!!HERE IS WHERE SCORE CUTOFF OF 50 IS SET!!!!!
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
                                    glycanMasses.Add(modMass);
                                }

                                if (modName.Contains("NGlycan"))
                                    Nlinked = true;

                                if (modName.Contains("OGlycan"))
                                    Olinked = true;

                                Modification modToAdd = new Modification(modMass, modName);
                                peptide.AddModification(modToAdd, modPosition);

                                if (modName.Contains("NGlycan"))
                                {
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


                        GlycoPSM psm = new GlycoPSM(peptide, peptideMonoMass, peptide.SequenceWithModifications, mods, glycans, glycanMasses, glycanPositionsList, uniprotID, PEP2D, logProb, score, deltaScore, deltaModScore, mzObs, charge, numberOfSites, ppmError,
                            obsMH, cleavage, proteinName, peptideStartPosition, scanTime, scanNumber, FDR2D, FDR2Dunique, qvalue2D, fragmentation, isGlycoPeptide, seenWithHCD, seenWithETD, NXSmotif, NXTmotif, isLocalized,
                            Nlinked, Olinked, matchedToUniprot, uniprotEvidenceType, uniprotEvidenceNumber);

                        if (isLocalized)
                        {
                            glycoPSMsLocalized.Add(psm);
                        }
                    }
                }
            }


            //look for fragments
            foreach (GlycoPSM glycoPSM in glycoPSMsLocalized)
            {
                if (glycoPSM.glycans.Count == 1)
                {

                    int glycoPosition = glycoPSM.glycanPositions[0];
                    double glycanMass = glycoPSM.glycanMasses[0];
                    string glycan = glycoPSM.glycans[0];

                    List<FragmentMatch> peptideFragments = new List<FragmentMatch>();
                    List<FragmentMatch> peptideNeutralLossFragments = new List<FragmentMatch>();
                    List<FragmentMatch> peptideFragmentsMustIncludeGlycan = new List<FragmentMatch>();
                    List<FragmentMatch> peptideIntactGlycans = new List<FragmentMatch>();
                    List<FragmentMatch> oxoniumIons = new List<FragmentMatch>();

                    int hexCount = 0;
                    int hexNAcCount = 0;
                    int fucCount = 0;
                    int neuAcCount = 0;
                    int neuGcCount = 0;

                    foreach (string glycanInList in glycoPSM.glycans)
                    {
                        hexCount += GetSugarCount(glycanInList, "Hex");
                        hexNAcCount += GetSugarCount(glycanInList, "HexNAc");
                        fucCount += GetSugarCount(glycanInList, "Fuc");
                        neuAcCount += GetSugarCount(glycanInList, "NeuAc");
                        neuGcCount += GetSugarCount(glycanInList, "NeuGc");
                    }


                    //generate theoretical fragments
                    List<Fragment> theoFragments = new List<Fragment>();

                    List<Fragment> bTypeFragments = glycoPSM.peptide.Fragment(FragmentTypes.b).ToList();
                    List<Fragment> yTypeFragments = glycoPSM.peptide.Fragment(FragmentTypes.y).ToList();
                    List<Fragment> cTypeFragments = glycoPSM.peptide.Fragment(FragmentTypes.c).ToList();
                    List<Fragment> zdotTypeFragments = glycoPSM.peptide.Fragment(FragmentTypes.zdot).ToList();

                    theoFragments.AddRange(bTypeFragments);
                    theoFragments.AddRange(yTypeFragments);

                    if (glycoPSM.fragmentation.Equals("ETD"))
                    {
                        theoFragments.AddRange(cTypeFragments);
                        theoFragments.AddRange(zdotTypeFragments);
                    }

                    bool notInterferingWithPrecursorPeaks = false;

                    //search for fragments in spectrum
                    foreach (Fragment theoFragment in theoFragments)
                    {
                        for (int i = glycoPSM.charge - 1; i > 0; i--)
                        {
                            double theoFragmentMZ = (theoFragment.MonoisotopicMass + (Constants.Hydrogen * i)) / ((double)i);

                            notInterferingWithPrecursorPeaks = NoPrecursorPeakInterference(theoFragmentMZ, glycoPSM.peptide, glycoPSM.charge, glycoPSM.fragmentation);

                            if (notInterferingWithPrecursorPeaks)
                            {
                                double theoFragmentIsoMZ = (theoFragmentMZ + (1 * (Constants.Hydrogen)) / ((double)i));
                                var range = DoubleRange.FromPPM(theoFragmentMZ, 20);

                                List<ThermoMzPeak> outPeaks = new List<ThermoMzPeak>();

                                if (rawFile.GetSpectrum(glycoPSM.scanNumber).TryGetPeaks(range, out outPeaks))
                                {
                                    var closestPeak = GetClosestPeak(outPeaks, theoFragmentMZ);
                                    double intensity = closestPeak.Intensity;

                                    var rangeForIsotope = DoubleRange.FromPPM(theoFragmentIsoMZ, 20);
                                    List<ThermoMzPeak> outPeaksIso = new List<ThermoMzPeak>();
                                    if (rawFile.GetSpectrum(glycoPSM.scanNumber).TryGetPeaks(rangeForIsotope, out outPeaksIso))
                                    {
                                        var closestPeakIso = GetClosestPeak(outPeaksIso, theoFragmentIsoMZ);
                                        intensity += closestPeakIso.Intensity;
                                    }

                                    double basePeak = rawFile.GetSpectrum(glycoPSM.scanNumber).GetBasePeakIntensity();

                                    if (closestPeak.Charge == i)
                                    {
                                        FragmentMatch fragmentMatch = new FragmentMatch(theoFragment.ToString(), theoFragment.Type.ToString(), theoFragment.Number, i, closestPeak.MZ, intensity);
                                        peptideFragments.Add(fragmentMatch);
                                        peptideFragmentsMustIncludeGlycan.Add(fragmentMatch);
                                    }
                                    else if (rawFile.GetSpectrum(glycoPSM.scanNumber).TryGetPeaks(rangeForIsotope, out outPeaksIso))
                                    {
                                        FragmentMatch fragmentMatch = new FragmentMatch(theoFragment.ToString(), theoFragment.Type.ToString(), theoFragment.Number, i, closestPeak.MZ, intensity);
                                        peptideFragments.Add(fragmentMatch);
                                        peptideFragmentsMustIncludeGlycan.Add(fragmentMatch);
                                    }
                                    else if (intensity >= (0.01 * basePeak) && closestPeak.Charge == 0)
                                    {
                                        FragmentMatch fragmentMatch = new FragmentMatch(theoFragment.ToString(), theoFragment.Type.ToString(), theoFragment.Number, i, closestPeak.MZ, intensity);
                                        peptideFragments.Add(fragmentMatch);
                                        peptideFragmentsMustIncludeGlycan.Add(fragmentMatch);
                                    }
                                }
                            }
                        }
                    }

                    //search for fragments without glycan that should have glycan
                    foreach(Fragment theoFragment in theoFragments)
                    {
                        LookForFragmentWithNoGlycan(rawFile, theoFragment, glycoPSM.charge, glycoPosition, glycanMass, peptideFragments, peptideNeutralLossFragments, glycoPSM.scanNumber, glycoPSM.peptide, glycoPSM.fragmentation);
                    }

                    //search for glycan peaks still attached to peptide
                    GetIntactPepGlycanFragments(glycoPSM, glycoPSM.charge, hexCount, hexNAcCount, fucCount, neuAcCount, neuGcCount, peptideIntactGlycans, rawFile, glycan, glycanMass);
                    GetOxoniumIons(glycoPSM, hexCount, hexNAcCount, fucCount, neuAcCount, neuGcCount, glycan, glycanMass, rawFile, oxoniumIons);

                    double sequenceCoverage = GetSequenceCoverage(glycoPSM, peptideFragments);
                    double sequenceCoverageOnlyGlycoFragments = GetSequenceCoverage(glycoPSM, peptideFragmentsMustIncludeGlycan);

                    bool anyFragmentsContainGlycan = SeeIfFragmentsHaveGlycan(glycoPSM, peptideFragmentsMustIncludeGlycan,"");
                    bool ntermFragmentsContainGlycan = SeeIfFragmentsHaveGlycan(glycoPSM, peptideFragmentsMustIncludeGlycan, "Nterm");
                    bool ctermFragmentsContainGlycan = SeeIfFragmentsHaveGlycan(glycoPSM, peptideFragmentsMustIncludeGlycan, "Cterm");

                    double glycanSequenceCoverage = GetGlycanSequenceCoverage(glycoPSM, peptideIntactGlycans, hexCount, hexNAcCount, fucCount, neuAcCount, neuGcCount);

                    string fragmentList = "";
                    double fragmentIntensity = 0;

                    double TIC = rawFile.GetSpectrum(glycoPSM.scanNumber).TotalIonCurrent;

                    double numberOfNeutralLossFragments = peptideNeutralLossFragments.Count;

                    outputEachFragment.Write(glycoPSM.peptide.ToString() + "\t" + glycoPSM.peptide.Sequence + "\t" + glycoPSM.glycans[0] + "\t" + fraction + "\t" + glycoPSM.scanNumber + "\t" + glycoPSM.fragmentation + "\t" +
                        glycoPSM.mzObs + "\t" + glycoPSM.charge + "\t" + glycoPSM.obsMH + "\t" + glycanMass + "\t" + glycoPSM.score + "\t" + glycoPSM.deltaModScore + "\t" + glycoPSM.FDR2D + "\t" + glycoPSM.PEP2D + "\t" +
                        sequenceCoverage + "\t" + sequenceCoverageOnlyGlycoFragments + "\t" + glycanSequenceCoverage + "\t" + peptideFragments.Count + "\t" + peptideNeutralLossFragments.Count + "\t" +
                        peptideIntactGlycans.Count + "\t" + oxoniumIons.Count + "\t");

                    outputPeptideFragments.Write(glycoPSM.peptide.ToString() + "\t" + glycoPSM.peptide.Sequence + "\t" + glycoPSM.glycans[0] + "\t" + fraction + "\t" + glycoPSM.scanNumber + "\t" + glycoPSM.fragmentation + "\t" +
                        glycoPSM.mzObs + "\t" + glycoPSM.charge + "\t" + glycoPSM.obsMH + "\t" + glycanMass + "\t" + glycoPSM.score + "\t" + glycoPSM.deltaModScore + "\t" + glycoPSM.FDR2D + "\t" + glycoPSM.PEP2D + "\t" +
                        sequenceCoverage + "\t" + sequenceCoverageOnlyGlycoFragments + "\t" + glycanSequenceCoverage + "\t" + peptideFragments.Count + "\t" + peptideNeutralLossFragments.Count + "\t");

                    outputGlycanFragments.Write(glycoPSM.peptide.ToString() + "\t" + glycoPSM.peptide.Sequence + "\t" + glycoPSM.glycans[0] + "\t" + fraction + "\t" + glycoPSM.scanNumber + "\t" + glycoPSM.fragmentation + "\t" +
                        glycoPSM.mzObs + "\t" + glycoPSM.charge + "\t" + glycoPSM.obsMH + "\t" + glycanMass + "\t" + glycoPSM.score + "\t" + glycoPSM.deltaModScore + "\t" + glycoPSM.FDR2D + "\t" + glycoPSM.PEP2D + "\t" +
                        sequenceCoverage + "\t" + sequenceCoverageOnlyGlycoFragments + "\t" + glycanSequenceCoverage + "\t" + peptideIntactGlycans.Count + "\t");

                    outputOxoniumFragments.Write(glycoPSM.peptide.ToString() + "\t" + glycoPSM.peptide.Sequence + "\t" + glycoPSM.glycans[0] + "\t" + fraction + "\t" + glycoPSM.scanNumber + "\t" + glycoPSM.fragmentation + "\t" +
                        glycoPSM.mzObs + "\t" + glycoPSM.charge + "\t" + glycoPSM.obsMH + "\t" + glycanMass + "\t" + glycoPSM.score + "\t" + glycoPSM.deltaModScore + "\t" + glycoPSM.FDR2D + "\t" + glycoPSM.PEP2D + "\t" +
                        sequenceCoverage + "\t" + sequenceCoverageOnlyGlycoFragments + "\t" + glycanSequenceCoverage + "\t" + oxoniumIons.Count + "\t");


                    foreach (FragmentMatch fragment in peptideFragments)
                    {
                        fragmentList += fragment.fragmentName + ";";
                        fragmentIntensity += fragment.fragmentSignal;
                        outputEachFragment.Write(fragment.fragmentName + "/" + fragment.fragmentNumber + "/" + fragment.fragmentMZ + "/" + fragment.fragmentCharge + "/" + fragment.fragmentSignal + "\t");
                        outputPeptideFragments.Write(fragment.fragmentName + "/" + fragment.fragmentNumber + "/" + fragment.fragmentMZ + "/" + fragment.fragmentCharge + "/" + fragment.fragmentSignal + "\t");
                    }

                    double fragmentsNoNeutralLossIntensity = 0;
                    foreach (FragmentMatch fragment in peptideFragmentsMustIncludeGlycan)
                    {
                        fragmentsNoNeutralLossIntensity += fragment.fragmentSignal;
                        //outputEachFragment.Write(fragment.fragmentName + "/" + fragment.fragmentNumber + "/" + fragment.fragmentMZ + "/" + fragment.fragmentCharge + "/" + fragment.fragmentSignal + "\t");
                    }

                    double glycanFragmentIntensity = 0;
                    string glycanFragmentList = "";
                    foreach (FragmentMatch fragment in peptideIntactGlycans)
                    {
                        glycanFragmentList += fragment.fragmentName + ";";
                        glycanFragmentIntensity += fragment.fragmentSignal;
                        outputEachFragment.Write(fragment.fragmentName + "/" + fragment.fragmentNumber + "/" + fragment.fragmentMZ + "/" + fragment.fragmentCharge + "/" + fragment.fragmentSignal + "\t");
                        outputGlycanFragments.Write(fragment.fragmentName + "/" + fragment.fragmentNumber + "/" + fragment.fragmentMZ + "/" + fragment.fragmentCharge + "/" + fragment.fragmentSignal + "\t");
                    }

                    double oxoniumIonIntensity = 0;
                    string oxoniumIonsList = "";
                    foreach (FragmentMatch fragment in oxoniumIons)
                    {
                        oxoniumIonsList += fragment.fragmentName + ";";
                        oxoniumIonIntensity += fragment.fragmentSignal;
                        outputEachFragment.Write(fragment.fragmentName + "/" + fragment.fragmentNumber + "/" + fragment.fragmentMZ + "/" + fragment.fragmentCharge + "/" + fragment.fragmentSignal + "\t");
                        outputOxoniumFragments.Write(fragment.fragmentName + "/" + fragment.fragmentNumber + "/" + fragment.fragmentMZ + "/" + fragment.fragmentCharge + "/" + fragment.fragmentSignal + "\t");
                    }

                    outputSummary.WriteLine(glycoPSM.peptide.ToString() + "\t" + glycoPSM.peptide.Sequence + "\t" + glycoPSM.glycans[0] + "\t" + fraction + "\t" + glycoPSM.scanNumber + "\t" + glycoPSM.fragmentation + "\t" +
                        glycoPSM.mzObs + "\t" + glycoPSM.charge + "\t" + glycoPSM.obsMH + "\t" + glycanMass + "\t" + glycoPSM.score + "\t" + glycoPSM.deltaModScore + "\t" + glycoPSM.FDR2D + "\t" + glycoPSM.PEP2D + "\t" +
                        sequenceCoverage + "\t" + sequenceCoverageOnlyGlycoFragments + "\t" + glycanSequenceCoverage + "\t" + peptideFragments.Count + "\t" + peptideNeutralLossFragments.Count + "\t" +
                        peptideIntactGlycans.Count + "\t" + oxoniumIons.Count + "\t" + anyFragmentsContainGlycan + "\t" + ntermFragmentsContainGlycan + "\t" + ctermFragmentsContainGlycan + "\t" + fragmentIntensity + "\t" +
                        fragmentsNoNeutralLossIntensity + "\t" + glycanFragmentIntensity + "\t" + oxoniumIonIntensity + "\t" + TIC + "\t" + fragmentList + "\t" + glycanFragmentList + "\t" + oxoniumIonsList);

                    outputEachFragment.WriteLine();
                    outputPeptideFragments.WriteLine();
                    outputGlycanFragments.WriteLine();
                    outputOxoniumFragments.WriteLine();
                }
            }

            outputSummary.Close();
            outputEachFragment.Close();
            outputPeptideFragments.Close();
            outputGlycanFragments.Close();
            outputOxoniumFragments.Close();
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

        public static int GetSugarCount(string glycan, string sugar)
        {
            int sugarCount = 0;

            string[] sugarsArray = glycan.Split(')');

            for (int i = 0; i < sugarsArray.Length; i++)
            {
                string[] glycanSep = sugarsArray[i].Split('(');
                if (glycanSep[0].Equals(sugar))
                {
                    sugarCount = Convert.ToInt32(glycanSep[1]);
                }
            }

            return sugarCount;
        }

        public static ThermoMzPeak GetClosestPeak(List<ThermoMzPeak> peaks, double theoMZ)
        {
            double diff = double.MaxValue;
            ThermoMzPeak returnPeak = null;
            foreach (var peak in peaks)
            {
                var currDiff = Math.Abs(peak.MZ - theoMZ);
                if (currDiff < diff)
                {
                    diff = currDiff;
                    returnPeak = peak;
                }
            }
            return returnPeak;
        }

        public static bool NoPrecursorPeakInterference(double mzOfFragment, Peptide peptide, int precursorChargeState, string fragmentation)
        {
            bool noInterference = false;

            bool noInterferenceIntactPrecursor = false;
            bool noInterferenceCRP = false;

            double intactPrecursorMZ = peptide.ToMz(precursorChargeState);

            if(mzOfFragment > (intactPrecursorMZ + 1) || mzOfFragment < (intactPrecursorMZ - 1))
            {
                noInterferenceIntactPrecursor = true;
            }

            if (fragmentation.Equals("ETD"))
            {
                for (int i = precursorChargeState - 1; i > 0; i--)
                {
                    double crp = (peptide.MonoisotopicMass + (Constants.Hydrogen*precursorChargeState))/ ((double) i);
                    if (mzOfFragment > (crp + 0.5) || mzOfFragment < (crp - 0.5))
                    {
                        noInterferenceCRP = true;
                    }
                    else
                    {
                        noInterferenceCRP = false;
                    }
                }
            }

            if(fragmentation.Equals("HCD"))
            {
                if(noInterferenceIntactPrecursor)
                {
                    noInterference = true;
                }
            }
            else if (fragmentation.Equals("ETD"))
            {
                if(noInterferenceIntactPrecursor && noInterferenceCRP)
                {
                    noInterference = true;
                }
            }

            return noInterference;
        }

        public static void LookForFragmentWithNoGlycan(ThermoRawFile rawFile, Fragment theoFragment, int charge, int glycanPosition, double glycanMass, List<FragmentMatch> peptideFragments, List<FragmentMatch> peptideNeutralLossFragments, int scanNumber, Peptide peptide, string fragmentation)
        {
            if (theoFragment.Type.Equals(FragmentTypes.b) || theoFragment.Type.Equals(FragmentTypes.c))
            {
                if (theoFragment.Number >= glycanPosition)
                {
                    for (int i = charge - 1; i > 0; i--)
                    {
                        double theoFragmentMZnoGlycan = (theoFragment.MonoisotopicMass + (Constants.Hydrogen * i) - glycanMass)  / ((double) i);

                        if (NoPrecursorPeakInterference(theoFragmentMZnoGlycan, peptide, charge, fragmentation))
                        {
                            double theoFragmentIsoMZnoGlycan = (theoFragmentMZnoGlycan + (1 * (Constants.Hydrogen)) / ((double)i));
                            var range = DoubleRange.FromPPM(theoFragmentMZnoGlycan, 20);

                            List<ThermoMzPeak> outPeaks = new List<ThermoMzPeak>();

                            if (rawFile.GetSpectrum(scanNumber).TryGetPeaks(range, out outPeaks))
                            {
                                var closestPeak = GetClosestPeak(outPeaks, theoFragmentMZnoGlycan);
                                double intensity = closestPeak.Intensity;

                                var rangeForIsotope = DoubleRange.FromPPM(theoFragmentIsoMZnoGlycan, 20);
                                List<ThermoMzPeak> outPeaksIso = new List<ThermoMzPeak>();
                                if (rawFile.GetSpectrum(scanNumber).TryGetPeaks(rangeForIsotope, out outPeaksIso))
                                {
                                    var closestPeakIso = GetClosestPeak(outPeaksIso, theoFragmentIsoMZnoGlycan);
                                    intensity += closestPeakIso.Intensity;
                                }

                                string fragmentName = "~" + theoFragment.ToString();

                                double basePeak = rawFile.GetSpectrum(scanNumber).GetBasePeakIntensity();

                                if (closestPeak.Charge == i)
                                {
                                    FragmentMatch fragmentMatch = new FragmentMatch(fragmentName, theoFragment.Type.ToString(), theoFragment.Number, i, closestPeak.MZ, intensity);
                                    peptideFragments.Add(fragmentMatch);
                                    peptideNeutralLossFragments.Add(fragmentMatch);
                                }
                                else if (rawFile.GetSpectrum(scanNumber).TryGetPeaks(rangeForIsotope, out outPeaksIso))
                                {
                                    FragmentMatch fragmentMatch = new FragmentMatch(fragmentName, theoFragment.Type.ToString(), theoFragment.Number, i, closestPeak.MZ, intensity);
                                    peptideFragments.Add(fragmentMatch);
                                    peptideNeutralLossFragments.Add(fragmentMatch);
                                }
                                else if (intensity >= (0.01 * basePeak) && closestPeak.Charge == 0)
                                {
                                    FragmentMatch fragmentMatch = new FragmentMatch(fragmentName, theoFragment.Type.ToString(), theoFragment.Number, i, closestPeak.MZ, intensity);
                                    peptideFragments.Add(fragmentMatch);
                                    peptideNeutralLossFragments.Add(fragmentMatch);
                                }
                            }
                        }

                        //Also look for 1 HexNAc still attached
                        ChemicalFormula hexNAc = new ChemicalFormula("C8O5NH13");

                        double theoFragmentMZplusHexNAc = (theoFragment.MonoisotopicMass + (Constants.Hydrogen * i) - glycanMass + hexNAc.MonoisotopicMass) / ((double)i);

                        if (NoPrecursorPeakInterference(theoFragmentMZplusHexNAc, peptide, charge, fragmentation))
                        {
                            double theoFragmentIsoMZplusHexNAc = (theoFragmentMZplusHexNAc + (1 * (Constants.Hydrogen)) / ((double)i));
                            var rangeHexNAcFragment = DoubleRange.FromPPM(theoFragmentMZplusHexNAc, 20);

                            List<ThermoMzPeak> outPeaksHexNAc = new List<ThermoMzPeak>();

                            if (rawFile.GetSpectrum(scanNumber).TryGetPeaks(rangeHexNAcFragment, out outPeaksHexNAc))
                            {
                                var closestPeak = GetClosestPeak(outPeaksHexNAc, theoFragmentMZplusHexNAc);
                                double intensity = closestPeak.Intensity;

                                var rangeForIsotopeHexNAc = DoubleRange.FromPPM(theoFragmentIsoMZplusHexNAc, 20);
                                List<ThermoMzPeak> outPeaksIsoHexNAc = new List<ThermoMzPeak>();
                                if (rawFile.GetSpectrum(scanNumber).TryGetPeaks(rangeForIsotopeHexNAc, out outPeaksIsoHexNAc))
                                {
                                    var closestPeakIso = GetClosestPeak(outPeaksIsoHexNAc, theoFragmentIsoMZplusHexNAc);
                                    intensity += closestPeakIso.Intensity;
                                }

                                string fragmentName = "~" + theoFragment.ToString() + "+203";

                                double basePeak = rawFile.GetSpectrum(scanNumber).GetBasePeakIntensity();

                                if (closestPeak.Charge == i)
                                {
                                    FragmentMatch fragmentMatch = new FragmentMatch(fragmentName, theoFragment.Type.ToString(), theoFragment.Number, i, closestPeak.MZ, intensity);
                                    peptideFragments.Add(fragmentMatch);
                                    peptideNeutralLossFragments.Add(fragmentMatch);
                                }
                                else if (rawFile.GetSpectrum(scanNumber).TryGetPeaks(rangeForIsotopeHexNAc, out outPeaksIsoHexNAc))
                                {
                                    FragmentMatch fragmentMatch = new FragmentMatch(fragmentName, theoFragment.Type.ToString(), theoFragment.Number, i, closestPeak.MZ, intensity);
                                    peptideFragments.Add(fragmentMatch);
                                    peptideNeutralLossFragments.Add(fragmentMatch);
                                }
                                else if (intensity >= (0.01 * basePeak) && closestPeak.Charge == 0)
                                {
                                    FragmentMatch fragmentMatch = new FragmentMatch(fragmentName, theoFragment.Type.ToString(), theoFragment.Number, i, closestPeak.MZ, intensity);
                                    peptideFragments.Add(fragmentMatch);
                                    peptideNeutralLossFragments.Add(fragmentMatch);
                                }
                            }
                        }
                    }
                }
            }
            else if (theoFragment.Type.Equals(FragmentTypes.y) || theoFragment.Type.Equals(FragmentTypes.zdot))
            {
                if (theoFragment.Number >= peptide.Length - glycanPosition)
                {
                    for (int i = charge - 1; i > 0; i--)
                    {
                        double theoFragmentMZnoGlycan = (theoFragment.MonoisotopicMass + (Constants.Hydrogen * i) - glycanMass) / ((double)i);

                        if (NoPrecursorPeakInterference(theoFragmentMZnoGlycan, peptide, charge, fragmentation))
                        {
                            double theoFragmentIsoMZnoGlycan = (theoFragmentMZnoGlycan + (1 * (Constants.Hydrogen)) / ((double)i));
                            var range = DoubleRange.FromPPM(theoFragmentMZnoGlycan, 20);

                            List<ThermoMzPeak> outPeaks = new List<ThermoMzPeak>();

                            if (rawFile.GetSpectrum(scanNumber).TryGetPeaks(range, out outPeaks))
                            {
                                var closestPeak = GetClosestPeak(outPeaks, theoFragmentMZnoGlycan);
                                double intensity = closestPeak.Intensity;

                                var rangeForIsotope = DoubleRange.FromPPM(theoFragmentIsoMZnoGlycan, 20);
                                List<ThermoMzPeak> outPeaksIso = new List<ThermoMzPeak>();
                                if (rawFile.GetSpectrum(scanNumber).TryGetPeaks(rangeForIsotope, out outPeaksIso))
                                {
                                    var closestPeakIso = GetClosestPeak(outPeaksIso, theoFragmentIsoMZnoGlycan);
                                    intensity += closestPeakIso.Intensity;
                                }

                                string fragmentName = "~" + theoFragment.ToString();

                                double basePeak = rawFile.GetSpectrum(scanNumber).GetBasePeakIntensity();

                                if (closestPeak.Charge == i)
                                {
                                    FragmentMatch fragmentMatch = new FragmentMatch(fragmentName, theoFragment.Type.ToString(), theoFragment.Number, i, closestPeak.MZ, intensity);
                                    peptideFragments.Add(fragmentMatch);
                                    peptideNeutralLossFragments.Add(fragmentMatch);
                                }
                                else if (rawFile.GetSpectrum(scanNumber).TryGetPeaks(rangeForIsotope, out outPeaksIso))
                                {
                                    FragmentMatch fragmentMatch = new FragmentMatch(fragmentName, theoFragment.Type.ToString(), theoFragment.Number, i, closestPeak.MZ, intensity);
                                    peptideFragments.Add(fragmentMatch);
                                    peptideNeutralLossFragments.Add(fragmentMatch);
                                }
                                else if (intensity >= (0.01 * basePeak) && closestPeak.Charge == 0)
                                {
                                    FragmentMatch fragmentMatch = new FragmentMatch(fragmentName, theoFragment.Type.ToString(), theoFragment.Number, i, closestPeak.MZ, intensity);
                                    peptideFragments.Add(fragmentMatch);
                                    peptideNeutralLossFragments.Add(fragmentMatch);
                                }
                            }
                        }

                        //Also look for 1 HexNAc still attached
                        ChemicalFormula hexNAc = new ChemicalFormula("C8O5NH13");

                        double theoFragmentMZplusHexNAc = (theoFragment.MonoisotopicMass + (Constants.Hydrogen * i) - glycanMass + hexNAc.MonoisotopicMass) / ((double)i);

                        if (NoPrecursorPeakInterference(theoFragmentMZplusHexNAc, peptide, charge, fragmentation))
                        {
                            double theoFragmentIsoMZplusHexNAc = (theoFragmentMZplusHexNAc + (1 * (Constants.Hydrogen)) / ((double)i));
                            var rangeHexNAcFragment = DoubleRange.FromPPM(theoFragmentMZplusHexNAc, 20);

                            List<ThermoMzPeak> outPeaksHexNAc = new List<ThermoMzPeak>();

                            if (rawFile.GetSpectrum(scanNumber).TryGetPeaks(rangeHexNAcFragment, out outPeaksHexNAc))
                            {
                                var closestPeak = GetClosestPeak(outPeaksHexNAc, theoFragmentMZplusHexNAc);
                                double intensity = closestPeak.Intensity;

                                var rangeForIsotopeHexNAc = DoubleRange.FromPPM(theoFragmentIsoMZplusHexNAc, 20);
                                List<ThermoMzPeak> outPeaksIsoHexNAc = new List<ThermoMzPeak>();
                                if (rawFile.GetSpectrum(scanNumber).TryGetPeaks(rangeForIsotopeHexNAc, out outPeaksIsoHexNAc))
                                {
                                    var closestPeakIso = GetClosestPeak(outPeaksIsoHexNAc, theoFragmentIsoMZplusHexNAc);
                                    intensity += closestPeakIso.Intensity;
                                }

                                string fragmentName = "~" + theoFragment.ToString() + "+203";

                                double basePeak = rawFile.GetSpectrum(scanNumber).GetBasePeakIntensity();

                                if (closestPeak.Charge == i)
                                {
                                    FragmentMatch fragmentMatch = new FragmentMatch(fragmentName, theoFragment.Type.ToString(), theoFragment.Number, i, closestPeak.MZ, intensity);
                                    peptideFragments.Add(fragmentMatch);
                                    peptideNeutralLossFragments.Add(fragmentMatch);
                                }
                                else if (rawFile.GetSpectrum(scanNumber).TryGetPeaks(rangeForIsotopeHexNAc, out outPeaksIsoHexNAc))
                                {
                                    FragmentMatch fragmentMatch = new FragmentMatch(fragmentName, theoFragment.Type.ToString(), theoFragment.Number, i, closestPeak.MZ, intensity);
                                    peptideFragments.Add(fragmentMatch);
                                    peptideNeutralLossFragments.Add(fragmentMatch);
                                }
                                else if (intensity >= (0.01 * basePeak) && closestPeak.Charge == 0)
                                {
                                    FragmentMatch fragmentMatch = new FragmentMatch(fragmentName, theoFragment.Type.ToString(), theoFragment.Number, i, closestPeak.MZ, intensity);
                                    peptideFragments.Add(fragmentMatch);
                                    peptideNeutralLossFragments.Add(fragmentMatch);
                                }
                            }
                        }

                    }
                }
            }
        }

        public static void GetIntactPepGlycanFragments(GlycoPSM glycoPSM, int charge, int hexCount, int hexNAcCount, int fucCount, int neuAcCount, int neuGcCount, List<FragmentMatch> peptideIntactGlycanFragments, ThermoRawFile rawFile, string glycanMod, double glycanMass)
        {
            Glycan glycan = new Glycan(glycanMod);
            glycan.AddHexNAc(hexNAcCount);
            glycan.AddHex(hexCount);
            glycan.AddFuc(fucCount);
            glycan.AddNeuAc(neuAcCount);
            glycan.AddNeuGc(neuGcCount);

            int totalSugarIntactGlycan = hexNAcCount + hexCount + fucCount + neuAcCount + neuGcCount;

            glycan.GenerateCombinations();

            var glycanPieces = glycan.AllFragments;

            var filteredPieces = GlycanMethods.GetValidGlycanStructures(glycanPieces);
            foreach (var piece in filteredPieces)
            {
                for (int i = charge; i > 0; i--)
                {
                    double peptidePlusGlycan_MZ = ((glycoPSM.peptide.MonoisotopicMass - glycanMass) + (Constants.Hydrogen * i) + GlycanMethods.GetGlycanMass(piece.Value)) / ((double)i);
                    if (peptidePlusGlycan_MZ < 2000)
                    {
                        double peptidePlusGlycan_MZiso1 = (peptidePlusGlycan_MZ + (1 * Constants.Hydrogen)) / ((double)i);
                        double peptidePlusGlycan_MZiso2 = (peptidePlusGlycan_MZ + (2 * Constants.Hydrogen)) / ((double)i);
                        double peptidePlusGlycan_MZiso3 = (peptidePlusGlycan_MZ + (3 * Constants.Hydrogen)) / ((double)i);

                        DoubleRange range = DoubleRange.FromPPM(peptidePlusGlycan_MZ, 40);
                        List<ThermoMzPeak> outPeaks = new List<ThermoMzPeak>();
                        if (rawFile.GetSpectrum(glycoPSM.scanNumber).TryGetPeaks(range, out outPeaks))
                        {
                            var closestPeak = GetClosestPeak(outPeaks, peptidePlusGlycan_MZ);
                            double intensity = closestPeak.Intensity;

                            DoubleRange rangeIso1 = DoubleRange.FromPPM(peptidePlusGlycan_MZiso1, 40);
                            List<ThermoMzPeak> outpeaksIso1 = new List<ThermoMzPeak>();
                            if (rawFile.GetSpectrum(glycoPSM.scanNumber).TryGetPeaks(rangeIso1, out outpeaksIso1))
                            {
                                var closestPeakIso1 = GetClosestPeak(outpeaksIso1, peptidePlusGlycan_MZiso1);
                                intensity += closestPeakIso1.Intensity;
                            }

                            DoubleRange rangeIso2 = DoubleRange.FromPPM(peptidePlusGlycan_MZiso2, 40);
                            List<ThermoMzPeak> outpeaksIso2 = new List<ThermoMzPeak>();
                            if (rawFile.GetSpectrum(glycoPSM.scanNumber).TryGetPeaks(rangeIso2, out outpeaksIso2))
                            {
                                var closestPeakIso2 = GetClosestPeak(outpeaksIso2, peptidePlusGlycan_MZiso2);
                                intensity += closestPeakIso2.Intensity;
                            }

                            DoubleRange rangeIso3 = DoubleRange.FromPPM(peptidePlusGlycan_MZiso3, 40);
                            List<ThermoMzPeak> outpeaksIso3 = new List<ThermoMzPeak>();
                            if (rawFile.GetSpectrum(glycoPSM.scanNumber).TryGetPeaks(rangeIso3, out outpeaksIso3))
                            {
                                var closestPeakIso3 = GetClosestPeak(outpeaksIso3, peptidePlusGlycan_MZiso3);
                                intensity += closestPeakIso3.Intensity;
                            }

                            double basePeak = rawFile.GetSpectrum(glycoPSM.scanNumber).GetBasePeakIntensity();

                            string fragmentName = "Pep+" + piece.Key + "_z" + i;

                            int sugarCountThisFragment = 0;
                            sugarCountThisFragment += GetSugarCount(piece.Key, "Hex");
                            sugarCountThisFragment += GetSugarCount(piece.Key, "HexNAc");
                            sugarCountThisFragment += GetSugarCount(piece.Key, "Fuc");
                            sugarCountThisFragment += GetSugarCount(piece.Key, "NeuAc");
                            sugarCountThisFragment += GetSugarCount(piece.Key, "NeuGc");

                            int sugarFragmentNumber = totalSugarIntactGlycan - sugarCountThisFragment;

                            if (closestPeak.Charge == i)
                            {
                                FragmentMatch fragmentMatch = new FragmentMatch(fragmentName, "Glyco", sugarFragmentNumber, i, peptidePlusGlycan_MZ, intensity);
                                peptideIntactGlycanFragments.Add(fragmentMatch);
                            }
                            else if(rawFile.GetSpectrum(glycoPSM.scanNumber).TryGetPeaks(rangeIso1, out outpeaksIso1))
                            {
                                FragmentMatch fragmentMatch = new FragmentMatch(fragmentName, "Glyco", sugarFragmentNumber, i, peptidePlusGlycan_MZ, intensity);
                                peptideIntactGlycanFragments.Add(fragmentMatch);
                            }
                            else if(intensity > (basePeak *0.01) && closestPeak.Charge == 0)
                            {
                                FragmentMatch fragmentMatch = new FragmentMatch(fragmentName, "Glyco", sugarFragmentNumber, i, peptidePlusGlycan_MZ, intensity);
                                peptideIntactGlycanFragments.Add(fragmentMatch);
                            }
                        }
                    }
                }
            }

            
        }

        public static void GetOxoniumIons(GlycoPSM glycoPSM, int hexCount, int hexNAcCount, int fucCount, int neuAcCount, int neuGcCount, string glycanMod, double glycanMass, ThermoRawFile rawFile, List<FragmentMatch> oxoniumIons)
        {
            Glycan glycan = new Glycan(glycanMod);
            glycan.AddHexNAc(hexNAcCount);
            glycan.AddHex(hexCount);
            glycan.AddFuc(fucCount);
            glycan.AddNeuAc(neuAcCount);
            glycan.AddNeuGc(neuGcCount);

            int totalSugarIntactGlycan = hexNAcCount + hexCount + fucCount + neuAcCount + neuGcCount;

            glycan.GenerateCombinations();

            var glycanPieces = glycan.AllFragments;

            var filteredPieces = GlycanMethods.GetValidGlycanStructures(glycanPieces);

            foreach (var piece in filteredPieces)
            {
                //Console.WriteLine(piece.Key + GlycanMethods.GetGlycanMZ(piece.Value, 1));
                FindOxoniumIon(GlycanMethods.GetGlycanMZ(piece.Value, 1), rawFile, oxoniumIons, glycoPSM.scanNumber, piece.Key);
            }

            FindOxoniumIon(366.14, rawFile, oxoniumIons, glycoPSM.scanNumber, "Hex(1)HexNAc(1)");
            FindOxoniumIon(528.19, rawFile, oxoniumIons, glycoPSM.scanNumber, "Hex(2)HexNAc(1)");
            FindOxoniumIon(690.245, rawFile, oxoniumIons, glycoPSM.scanNumber, "Hex(3)HexNAc(1)");

            if (hexCount > 0)
            {
                FindOxoniumIon(163.06, rawFile, oxoniumIons, glycoPSM.scanNumber, "Hex(1)");
            }
            if (fucCount > 0)
            {
                FindOxoniumIon(147.065, rawFile, oxoniumIons, glycoPSM.scanNumber, "Fuc(1)");
            }
            if (neuAcCount > 0)
            {
                FindOxoniumIon(274.092, rawFile, oxoniumIons, glycoPSM.scanNumber, "NeuAcWaterLoss");
                FindOxoniumIon(292.103, rawFile, oxoniumIons, glycoPSM.scanNumber, "NeuAc(1)");
                FindOxoniumIon(495.1382, rawFile, oxoniumIons, glycoPSM.scanNumber, "HexNAc(1)NeuAc(1)");
                FindOxoniumIon(657.23, rawFile, oxoniumIons, glycoPSM.scanNumber, "Hex(1)HexNAc(1)NeuAc(1)");

            }
            if (neuGcCount > 0)
            {
                FindOxoniumIon(290.087, rawFile, oxoniumIons, glycoPSM.scanNumber, "NeuGcWaterLoss");
                FindOxoniumIon(308.098, rawFile, oxoniumIons, glycoPSM.scanNumber, "NeuGc(1)");
            }
            
        }

        public static void FindOxoniumIon(double mz, ThermoRawFile rawFile, List<FragmentMatch> list, int scanNumber, string name)
        {
            double mzIso = mz + Constants.Hydrogen;

            DoubleRange mzRange = DoubleRange.FromPPM(mz, 20);
            List<ThermoMzPeak> outPeaks = new List<ThermoMzPeak>();

            if (rawFile.GetSpectrum(scanNumber).TryGetPeaks(mzRange, out outPeaks))
            {
                var closestPeak = GetClosestPeak(outPeaks, mz);
                double intensity = closestPeak.Intensity;

                DoubleRange isoMZrange = DoubleRange.FromPPM(mzIso, 20);
                List<ThermoMzPeak> outPeaksIso = new List<ThermoMzPeak>();
                if (rawFile.GetSpectrum(scanNumber).TryGetPeaks(isoMZrange, out outPeaksIso))
                {
                    var closestPeakIso = GetClosestPeak(outPeaksIso, mzIso);
                    intensity += closestPeakIso.Intensity;
                }

                double basePeak = rawFile.GetSpectrum(scanNumber).GetBasePeakIntensity();

                if(intensity > (0.01 * basePeak))
                {
                    FragmentMatch match = new FragmentMatch(name, "Oxonium", 0, 1, mz, intensity);
                    list.Add(match);
                }
            }

        }

        public static double GetSequenceCoverage(GlycoPSM glycoPSM, List<FragmentMatch> fragments)
        {
            double sequenceCoverage = 0;

            char[] peptidePostitionArrary = new char[glycoPSM.peptide.Length-1];

            foreach (FragmentMatch fragment in fragments)
            {
                if (fragment.fragmentType.Equals("b") || fragment.fragmentType.Equals("c"))
                {
                    int bondNumberBroken = fragment.fragmentNumber;
                    int positionInArray = bondNumberBroken - 1;
                    peptidePostitionArrary[positionInArray] = 'X';
                }
                if (fragment.fragmentType.Equals("y") || fragment.fragmentType.Equals("zdot"))
                {
                    int bondNumberBroken = glycoPSM.peptide.Length - fragment.fragmentNumber;
                    int positionInArray = bondNumberBroken - 1;
                    peptidePostitionArrary[positionInArray] = 'X';
                }
            }

            int numberOfBondsBroken = 0;
            for(int i = 0; i < peptidePostitionArrary.Length; i ++)
            {
                if (peptidePostitionArrary[i].Equals('X'))
                    numberOfBondsBroken++;
            }
            sequenceCoverage = (double)numberOfBondsBroken / (double) peptidePostitionArrary.Length;

            return sequenceCoverage;
        }

        public static double GetGlycanSequenceCoverage(GlycoPSM glycoPSM, List<FragmentMatch> fragments, int hexCount, int hexNAcCount, int fucCount, int neuAcCount, int neuGcCount)
        {
            double sequenceCoverage = 0;
            Glycan glycan = new Glycan(glycoPSM.glycans[0]);
            glycan.AddHexNAc(hexNAcCount);
            glycan.AddHex(hexCount);
            glycan.AddFuc(fucCount);
            glycan.AddNeuAc(neuAcCount);
            glycan.AddNeuGc(neuGcCount);

            glycan.GenerateCombinations();

            var glycanPieces = glycan.AllFragments;

            var filteredPieces = GlycanMethods.GetValidGlycanStructures(glycanPieces);

            HashSet<string> theorecticalFragments = new HashSet<string>();

            //Console.WriteLine(glycoPSM.glycans[0]);
            foreach (var piece in filteredPieces)
            {
                string fragment = "Pep+" + piece.Key;
                theorecticalFragments.Add(fragment);
            }

            HashSet<string> matches = new HashSet<string>();

            foreach (FragmentMatch match in fragments)
            {
                string[] glycanFrag = match.fragmentName.Split('_');
                matches.Add(glycanFrag[0]);
            }

            int numberOfGlycoFragmentsExplained = 0;
            int numberOfGlycoFragmentsTotal = theorecticalFragments.Count;

            foreach (string theoFrag in theorecticalFragments)
            {
                if (matches.Contains(theoFrag))
                {
                    numberOfGlycoFragmentsExplained++;
                }
            }

            sequenceCoverage = (double)numberOfGlycoFragmentsExplained / (double)numberOfGlycoFragmentsTotal;

            return sequenceCoverage;
        }

        public static bool SeeIfFragmentsHaveGlycan(GlycoPSM glycoPSM, List<FragmentMatch> fragments, string terminus)
        {
            bool returnBoolAll = false;

            bool returnBoolNterm = false;
            bool returnBoolCterm = false;



            int glycanPosition = glycoPSM.glycanPositions[0];
            //Console.WriteLine(glycoPSM.peptide);
            //Console.WriteLine(glycanPosition);

            foreach (FragmentMatch fragment in fragments)
            {
                //Console.WriteLine(fragment.fragmentName);
                if (fragment.fragmentType.Equals("b") || fragment.fragmentType.Equals("c"))
                {
                    if (fragment.fragmentNumber >= glycanPosition)
                    {
                        returnBoolAll = true;
                        returnBoolNterm = true;
                    }
                }

                if (fragment.fragmentType.Equals("y") || fragment.fragmentType.Equals("zdot"))
                {
                    if (fragment.fragmentNumber > (glycoPSM.peptide.Length - glycanPosition))
                    {
                        returnBoolAll = true;
                        returnBoolCterm = true;
                    }
                }
            }

            //Console.WriteLine(returnBoolAll);
            //Console.WriteLine(returnBoolNterm);
            //Console.WriteLine(returnBoolCterm);
            //Console.ReadKey();
            if (terminus.Equals("NTerm"))
            {
                return returnBoolNterm;
            }
            else if (terminus.Equals("Cterm"))
            {
                return returnBoolCterm;
            }
            else
            {
                return returnBoolAll;
            }
        }


        /**
        public static void GetIntactPepGlycanFragments(GlycoPSM glycoPSM, int charge, int hexCount, int hexNAcCount, int fucCount, int neuAcCount, int neuGcCount, List<FragmentMatch> peptideIntactGlycanFragments, bool highMannoseGlycan, ThermoRawFile rawFile, string glycan)
        {
            ChemicalFormula hexNAc = new ChemicalFormula("C8O5NH13");
            ChemicalFormula hex = new ChemicalFormula("C6O5H10");
            ChemicalFormula fuc = new ChemicalFormula("C6O4H10");
            ChemicalFormula neuAc = new ChemicalFormula("C11O8NH17");
            ChemicalFormula neuGc = new ChemicalFormula("C11O9NH17");

            bool hexNAcPeakFound = false;
            bool hexPeakFound = false;

            for(int i = charge - 1; i > 0; i--)
            {
                double peptidePlusHexNAc1_MZ = (glycoPSM.peptide.MonoisotopicMass + (Constants.Hydrogen * i) + hexNAc.MonoisotopicMass) / ((double)i);
                double peptidePlusHexNAc1_MZ_iso = (peptidePlusHexNAc1_MZ + (1 * Constants.Hydrogen))/ ((double) i);

                string fragmentName = "Pep+HexNAc(1)";

                hexNAcPeakFound = FindFragments(peptidePlusHexNAc1_MZ, peptidePlusHexNAc1_MZ_iso, rawFile, glycoPSM.scanNumber, fragmentName, charge, peptideIntactGlycanFragments);
            }

            if (hexNAcPeakFound)
            {
                hexNAcCount -= 1;
                hexNAcPeakFound = false;
            }

            if (hexNAcCount > 0)
            {
                for (int i = charge - 1; i > 0; i--)
                {
                    double peptidePlusHexNAc2_MZ = (glycoPSM.peptide.MonoisotopicMass + (Constants.Hydrogen * i) + (2 * hexNAc.MonoisotopicMass)) / ((double)i);
                    double peptidePlusHexNAc2_MZ_iso = (peptidePlusHexNAc2_MZ + (1 * Constants.Hydrogen)) / ((double)i);

                    string fragmentName = "Pep+HexNAc(2)";

                    hexNAcPeakFound = FindFragments(peptidePlusHexNAc2_MZ, peptidePlusHexNAc2_MZ_iso, rawFile, glycoPSM.scanNumber, fragmentName, charge, peptideIntactGlycanFragments);
                }

                if (hexNAcPeakFound)
                {
                    hexNAcCount -= 1;
                    hexNAcPeakFound = false;
                }
            }

            if (hexCount > 0)
            {
                for (int i = charge - 1; i > 0; i--)
                {
                    double peptidePlusHexNAc2Hex1_MZ = (glycoPSM.peptide.MonoisotopicMass + (Constants.Hydrogen * i) + (2 * hexNAc.MonoisotopicMass) + (1* hex.MonoisotopicMass)) / ((double)i);
                    double peptidePlusHexNAc2Hex1_MZ_iso = (peptidePlusHexNAc2Hex1_MZ + (1 * Constants.Hydrogen)) / ((double)i);

                    string fragmentName = "Pep+HexNAc(2)Hex(1)";

                    hexPeakFound = FindFragments(peptidePlusHexNAc2Hex1_MZ, peptidePlusHexNAc2Hex1_MZ_iso, rawFile, glycoPSM.scanNumber, fragmentName, charge, peptideIntactGlycanFragments);
                }

                if (hexPeakFound)
                {
                    hexCount -= 1;
                    hexPeakFound = false;
                }
            }

            if (hexCount > 0)
            {
                for (int i = charge - 1; i > 0; i--)
                {
                    double peptidePlusHexNAc2Hex2_MZ = (glycoPSM.peptide.MonoisotopicMass + (Constants.Hydrogen * i) + (2 * hexNAc.MonoisotopicMass) + (2 * hex.MonoisotopicMass)) / ((double)i);
                    double peptidePlusHexNAc2Hex2_MZ_iso = (peptidePlusHexNAc2Hex2_MZ + (1 * Constants.Hydrogen)) / ((double)i);

                    string fragmentName = "Pep+HexNAc(2)Hex(2)";

                    hexPeakFound = FindFragments(peptidePlusHexNAc2Hex2_MZ, peptidePlusHexNAc2Hex2_MZ_iso, rawFile, glycoPSM.scanNumber, fragmentName, charge, peptideIntactGlycanFragments);
                }

                if (hexPeakFound)
                {
                    hexCount -= 1;
                    hexPeakFound = false;
                }
            }

            if (hexCount > 0)
            {
                for (int i = charge - 1; i > 0; i--)
                {
                    double peptidePlusHexNAc2Hex3_MZ = (glycoPSM.peptide.MonoisotopicMass + (Constants.Hydrogen * i) + (2 * hexNAc.MonoisotopicMass) + (3 * hex.MonoisotopicMass)) / ((double)i);
                    double peptidePlusHexNAc2Hex3_MZ_iso = (peptidePlusHexNAc2Hex3_MZ + (1 * Constants.Hydrogen)) / ((double)i);

                    string fragmentName = "Pep+HexNAc(2)Hex(3)";

                    hexPeakFound = FindFragments(peptidePlusHexNAc2Hex3_MZ, peptidePlusHexNAc2Hex3_MZ_iso, rawFile, glycoPSM.scanNumber, fragmentName, charge, peptideIntactGlycanFragments);
                }

                if (hexPeakFound)
                {
                    hexCount -= 1;
                    hexPeakFound = false;
                }
            }

            if (highMannoseGlycan)
            {
                int totalMannoseCount = hexCount;
                for (int j = 4; j <= hexCount; j++)
                {
                    for (int i = charge - 1; i > 0; i--)
                    {
                        double peptidePlusHexNAc2HexX_MZ = (glycoPSM.peptide.MonoisotopicMass + (Constants.Hydrogen * i) + (2 * hexNAc.MonoisotopicMass) + (j * hex.MonoisotopicMass)) / ((double)i);
                        double peptidePlusHexNAc2HexX_MZ_iso = (peptidePlusHexNAc2HexX_MZ + (1 * Constants.Hydrogen)) / ((double)i);

                        string fragmentName = "Pep+HexNAc(2)Hex(" + j + ")";

                        hexPeakFound = FindFragments(peptidePlusHexNAc2HexX_MZ, peptidePlusHexNAc2HexX_MZ_iso, rawFile, glycoPSM.scanNumber, fragmentName, charge, peptideIntactGlycanFragments);

                        //Console.WriteLine(glycoPSM.glycans[0]);
                        if(hexPeakFound)
                        {
                            Console.WriteLine(fragmentName);
                            Console.ReadKey();
                        }
                    }
                }
            }
            

            
        }
         * */

        /*
        public static bool FindFragments(double MZ, double isoMZ, ThermoRawFile rawFile, int specNumber, string name, int charge, List<FragmentMatch> list)
        {
            bool found = false;
            var range = DoubleRange.FromPPM(MZ, 20);
            List<ThermoMzPeak> outpeaks = new List<ThermoMzPeak>();
            if(rawFile.GetSpectrum(specNumber).TryGetPeaks(range, out outpeaks))
            {
                var closestPeak = GetClosestPeak(outpeaks, MZ);
                double intensity = closestPeak.Intensity;

                var rangeIso = DoubleRange.FromPPM(isoMZ, 20);
                List<ThermoMzPeak> outPeaksIso = new List<ThermoMzPeak>();

                if (rawFile.GetSpectrum(specNumber).TryGetPeaks(rangeIso, out outPeaksIso))
                {
                    var closestPeakIso = GetClosestPeak(outPeaksIso, isoMZ);
                    intensity += closestPeakIso.Intensity;
                }

                FragmentMatch fragmentFound = new FragmentMatch(name, "glyco", 0, charge, closestPeak.MZ, intensity);
                list.Add(fragmentFound);
                found = true;
            }
            return found;
        }
         * */
    }
}
