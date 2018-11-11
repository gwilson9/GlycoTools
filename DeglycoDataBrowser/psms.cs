using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DeglycoDataBrowser
{
    class psms
    {
        public string spectrumNumber;
        public string rawFileName;
        public string mass;
        public string protein;
        public string sequence;
        public string charge;
        public List<int> deglycoMods; //Position within protein sequence


        public psms(string specNum, string mass, string prot, string seq, string chrg, string mod, string start, string rawFileName)
        {
            this.spectrumNumber = specNum;
            this.mass = mass;
            this.protein = prot;
            this.sequence = seq;
            this.charge = chrg;
            this.rawFileName = rawFileName.Split('.')[0] + ".raw";
            this.deglycoMods = new List<int>();

            List<string> mods = mod.Split(',').ToList();
            foreach (string modd in mods)
            {
                if (!String.IsNullOrEmpty(modd))
                {
                    int location = Int32.Parse(modd.Split(':')[1]);

                    if (modd.Contains("deamid") && seq[location - 1].Equals('n'))
                    {
                        deglycoMods.Add(Int32.Parse(start) + location - 1);
                    }
                }
            }
        }

        public bool Equals(psms psm1, psms psm2)
        {
            if (psm1.sequence.Equals(psm2.sequence) && psm1.protein.Equals(psm2.protein))
                return true;

            return false;
        }

        public int GetHashCode(psms psm)
        {
            return psm.sequence.GetHashCode() + psm.protein.GetHashCode();
        }
    }
}
