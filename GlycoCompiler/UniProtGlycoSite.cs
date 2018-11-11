using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GlycoCompiler
{
    class UniProtGlycoSite
    {
        public string uniprotID { get; set; }
        public string defLine { get; set; }
        public double proteinMW { get; set; }
        public int glycoSite { get; set; }
        public string glycoType { get; set; }
        public string evidenceType { get; set; }
        public int evidenceNumber { get; set; }

        public UniProtGlycoSite(string uniprotID, string defLine, double proteinMw, int glycoSite, string glycoType, string evidenceType, int evidenceNumber)
        {
            this.uniprotID = uniprotID;
            this.defLine = defLine;
            this.proteinMW = proteinMW;
            this.glycoSite = glycoSite;
            this.glycoType = glycoType;
            this.evidenceType = evidenceType;
            this.evidenceNumber = evidenceNumber;
        }
    }
}
