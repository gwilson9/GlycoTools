using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ScanAssigner
{
    class HighlightEventArgs : EventArgs
    {
        public string rawFile { get; private set; }

        public string resultsFile { get; private set; }

        public HighlightEventArgs(string resultsFile, string rawFile)
        {
            this.rawFile = rawFile;
            this.resultsFile = resultsFile;
        }
    }
}
