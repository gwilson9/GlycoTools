using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using CSMSL.Proteomics;

namespace ScanAssigner
{
    class GlycoPSM
    {
        public Peptide peptide { get; set; }
        public double peptideMonoMass { get; set; }
        public string sequenceWithMods { get; set; }
        public List<string> mods { get; set; }
        public List<string> glycans { get; set; }
        public List<double> glycanMasses { get; set; }
        public List<int> glycanPositions { get; set; }
        public string uniprotID { get; set; }
        public double PEP2D { get; set; }
        public double logProb { get; set; }
        public double score { get; set; }
        public double deltaScore { get; set; }
        public double deltaModScore { get; set; }
        public double mzObs { get; set; }
        public int charge { get; set; }
        public int numberOfSites { get; set; }
        public double ppmError { get; set; }
        public double obsMH { get; set; }
        public string cleavage { get; set; }
        public string proteinName { get; set; }
        public int peptideStartPosition { get; set; }
        public double scanTime { get; set; }
        public int scanNumber { get; set; }
        public double FDR2D { get; set; }
        public double FDR2Dunique { get; set; }
        public double qvalue2D { get; set; }
        public string fragmentation { get; set; }
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

        public GlycoPSM(Peptide peptide, double peptideMonoMass, string sequenceWithMods, List<string> mods, List<string> glycans, List<double> glycanMasses, List<int> glycanPositions, string uniprotID, double PEP2D, double logProb, double score,
            double deltaScore, double deltaModScore, double mzObs, int charge, int numberOfSites, double ppmError, double obsMH, string cleavage, string proteinName,
            int peptideStartPosition, double scanTime, int scanNumber, double FDR2D, double FDR2Dunique, double qvalue2D, string fragmentation, bool isGlycoPeptide,
            bool seenWithHCD, bool seenWithETD, bool NXSmotif, bool NXTmotif, bool isLocalized, bool Nlinked, bool Olinked, bool matchedToUniprot, string evidenceType, int evidenceNumber)
        {
            this.peptide = peptide;
            this.peptideMonoMass = peptideMonoMass;
            this.sequenceWithMods = sequenceWithMods;
            this.mods = mods;
            this.glycans = glycans;
            this.glycanMasses = glycanMasses;
            this.glycanPositions = glycanPositions;
            this.uniprotID = uniprotID;
            this.PEP2D = PEP2D;
            this.logProb = logProb;
            this.score = score;
            this.deltaScore = deltaScore;
            this.deltaModScore = deltaModScore;
            this.mzObs = mzObs;
            this.charge = charge;
            this.numberOfSites = numberOfSites;
            this.ppmError = ppmError;
            this.obsMH = obsMH;
            this.cleavage = cleavage;
            this.proteinName = proteinName;
            this.peptideStartPosition = peptideStartPosition;
            this.scanTime = scanTime;
            this.scanNumber = scanNumber;
            this.FDR2D = FDR2D;
            this.FDR2Dunique = FDR2Dunique;
            this.qvalue2D = qvalue2D;
            this.fragmentation = fragmentation;
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
        }
    }
}
