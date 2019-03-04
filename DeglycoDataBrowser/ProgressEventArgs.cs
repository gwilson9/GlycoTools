using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DeglycoDataBrowser
{
    class ProgressEventArgs : EventArgs
    {
        public double Progress { get; private set; }

        public ProgressEventArgs(double progress)
        {
            Progress = progress;
        }
    }
}
