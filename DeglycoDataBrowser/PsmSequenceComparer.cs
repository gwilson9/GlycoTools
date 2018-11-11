using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DeglycoDataBrowser
{
    class PsmSequenceComparer : IEqualityComparer<psms>
    {
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
