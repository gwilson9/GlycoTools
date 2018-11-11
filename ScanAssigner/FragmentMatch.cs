using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
namespace ScanAssigner
{
    class FragmentMatch
    {
        public string fragmentName { get; set; }
        public string fragmentType { get; set; }
        public int fragmentNumber { get; set; }
        public int fragmentCharge { get; set; }
        public double fragmentMZ { get; set; }
        public double fragmentSignal { get; set; }

        public FragmentMatch(string fragmentName, string fragmentType, int fragmentNumber, int fragmentCharge, double fragmentMZ, double fragmentSignal)
        {
            this.fragmentName = fragmentName;
            this.fragmentType = fragmentType;
            this.fragmentNumber = fragmentNumber;
            this.fragmentCharge = fragmentCharge;
            this.fragmentMZ = fragmentMZ;
            this.fragmentSignal = fragmentSignal;
        }
    }
}
