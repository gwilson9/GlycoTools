using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using CSMSL.Proteomics;

namespace GlycoCompiler
{
    class GlycoPSM
    {
        public Peptide peptide { get; set; }
        public double peptideMonoMass { get; set; }
        public string sequenceWithMods { get; set; }
        public List<string> mods { get; set; }
        public List<string> glycans { get; set; }
        public List<int> glycanPositions { get; set; }
        public string uniprotID { get; set; }
        public double PEP2D { get; set; }
        public double PEP1D { get; set; }
        public double logProb { get; set; }
        public double score { get; set; }
        public double deltaScore { get; set; }
        public double deltaModScore { get; set; }
        public double mzObs { get; set; }
        public double mzCalc { get; set; }
        public int charge { get; set; }
        public int numberOfSites { get; set; }
        public double ppmError { get; set; }
        public double obsMH { get; set; }
        public double calcMH { get; set; }
        public string cleavage { get; set; }
        public string proteinName { get; set; }
        public int peptideStartPosition { get; set; }
        public double scanTime { get; set; }
        public int scanNumber { get; set; }
        public double FDR2D { get; set; }
        public double FDR1D { get; set; }
        public double FDR2Dunique { get; set; }
        public double FDR1Dunique { get; set; }
        public double qvalue2D { get; set; }
        public double qvalue1D { get; set; }
        public string fragmentation { get; set; }
        public int masterScan { get; set; }
        public double peak126 { get; set; }
        public double peak138 { get; set; }
        public double peak144 { get; set; }
        public double peak168 { get; set; }
        public double peak186 { get; set; }
        public double peak204 { get; set; }
        public double peak274 { get; set; }
        public double peak292 { get; set; }
        public double peak366 { get; set; }
        public double GlcNAcGalNAcratio { get; set; }
        public bool isGlycoPeptide { get; set; }
        public bool seenWithHCD { get; set; }
        public bool seenWithETD { get; set; }
        public bool NXSmotif { get; set; }
        public bool NXTmotif { get; set; }
        public bool isLocalized { get; set; }
        public bool Nlinked { get; set; }
        public bool Olinked { get; set; }
        public bool matchedToUniprot { get; set; }
        public string evidenceType { get; set; }
        public int evidenceNumber { get; set; }
        public double intensity;
        public string byonicIntensity;
        public List<string> glycanTypes;
        public string file;

        public GlycoPSM(Peptide peptide, double peptideMonoMass, string sequenceWithMods, List<string> mods, 
            List<string> glycans, List<int> glycanPositions, string uniprotID, double PEP2D, double PEP1D, double logProb, double score,
            double deltaScore, double deltaModScore, double mzObs, double mzCalc, int charge, int numberOfSites, double ppmError, 
            double obsMH, double calcMH, string cleavage, string proteinName, int peptideStartPosition, double scanTime, 
            int scanNumber, double FDR2D, double FDR1D, double FDR2Dunique, double FDR1Dunique, double qvalue2D, double qvalue1D, string fragmentation,
            int masterScan, double peak126, double peak138, double peak144, double peak168, double peak186, double peak204, double peak274, 
            double peak292, double peak366, double GlcNAcGalNAcratio, bool isGlycoPeptide,bool seenWithHCD, bool seenWithETD, bool NXSmotif, 
            bool NXTmotif, bool isLocalized, bool Nlinked, bool Olinked, bool matchedToUniprot, string evidenceType, int evidenceNumber, 
            double intensity, string byonicIntensity, string file)
        {
            this.peptide = peptide;
            this.peptideMonoMass = peptideMonoMass;
            this.sequenceWithMods = sequenceWithMods;
            this.mods = mods;
            this.glycans = glycans;
            this.glycanPositions = glycanPositions;
            this.uniprotID = uniprotID;
            this.PEP2D = PEP2D;
            this.PEP1D = PEP1D;
            this.logProb = logProb;
            this.score = score;
            this.deltaScore = deltaScore;
            this.deltaModScore = deltaModScore;
            this.mzObs = mzObs;
            this.mzCalc = mzCalc;
            this.charge = charge;
            this.numberOfSites = numberOfSites;
            this.ppmError = ppmError;
            this.obsMH = obsMH;
            this.calcMH = calcMH;
            this.cleavage = cleavage;
            this.proteinName = proteinName;
            this.peptideStartPosition = peptideStartPosition;
            this.scanTime = scanTime;
            this.scanNumber = scanNumber;
            this.FDR2D = FDR2D;
            this.FDR1D = FDR1D;
            this.FDR2Dunique = FDR2Dunique;
            this.FDR1Dunique = FDR1Dunique;
            this.qvalue2D = qvalue2D;
            this.qvalue1D = qvalue1D;
            this.fragmentation = fragmentation;
            this.masterScan = masterScan;
            this.peak126 = peak126;
            this.peak138 = peak138;
            this.peak144 = peak144;
            this.peak168 = peak168;
            this.peak186 = peak186;
            this.peak204 = peak204;
            this.peak274 = peak274;
            this.peak292 = peak292;
            this.peak366 = peak366;
            this.GlcNAcGalNAcratio = GlcNAcGalNAcratio;
            this.isGlycoPeptide = isGlycoPeptide;
            this.seenWithHCD = seenWithHCD;
            this.seenWithETD = seenWithETD;
            this.NXSmotif = NXSmotif;
            this.NXTmotif = NXTmotif;
            this.isLocalized = isLocalized;
            this.Nlinked = Nlinked;
            this.Olinked = Olinked;
            this.matchedToUniprot = matchedToUniprot;
            this.evidenceType = evidenceType;
            this.evidenceNumber = evidenceNumber;
            this.intensity = intensity;
            this.byonicIntensity = byonicIntensity;
            this.glycanTypes = new List<string> ();
            this.file = file;

            foreach(string glycan in glycans)
            {
                string type = "";

                for (int i = 4; i < 13; i++)
                {
                    string name = "HexNAc(2)Hex(" + i + ")";

                    if (glycan.Equals(name))
                    {
                        type = "Mann";
                        //typeNumber = 2;
                    }
                }

                if (glycan.Equals("HexNAc(1)") || glycan.Equals("HexNAc(1)Fuc(1)") || glycan.Equals("HexNAc(2)") || glycan.Equals("HexNAc(2)Fuc(1)") ||
                    glycan.Equals("HexNAc(2)Hex(1)") || glycan.Equals("HexNAc(2)Hex(1)Fuc(1)") || glycan.Equals("HexNAc(2)Hex(2)") ||
                    glycan.Equals("HexNAc(2)Hex(2)Fuc(1)") || glycan.Equals("HexNAc(2)Hex(3)") || glycan.Equals("HexNAc(2)Hex(3)Fuc(1)"))
                {
                    type = "Trunc";
                    //typeNumber = 1;
                }

                if (glycan.Contains("NeuAc"))
                {
                    type = "NeuAc";
                    //typeNumber = 6;
                }

                if (String.IsNullOrEmpty(type) && glycan.Contains("Fuc") && !glycan.Contains("NeuAc"))
                {
                    type = "Fuc";
                    //typeNumber = 5;
                }

                if (glycan.Equals("HexNAc(2)Hex(6)Phospho(1)"))
                {
                    type = "MannPhospho";
                    //typeNumber = 3;
                }

                if (String.IsNullOrEmpty(type))
                {
                    type = "Other";
                    //typeNumber = 4;
                }

                glycanTypes.Add(type);

            }
        }
    }
}
