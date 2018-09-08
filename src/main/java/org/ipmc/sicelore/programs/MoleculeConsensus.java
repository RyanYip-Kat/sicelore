package org.ipmc.sicelore.programs;

import java.io.*;
import java.util.*;
import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import org.ipmc.sicelore.utils.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "Procude molecule consensus sequence", oneLineSummary = "Procude molecule consensus sequence", programGroup = org.ipmc.sicelore.cmdline.SiCeLoRe.class)
@DocumentedFeature
public class MoleculeConsensus extends CommandLineProgram {

    private final Log log;
    private htsjdk.samtools.util.ProgressLogger pl;

    @Argument(shortName = "I", doc = "The input SAM or BAM file to analyze")
    public File INPUT;
    @Argument(shortName = "O", doc = "The output .fasta file of consensus sequences")
    public File OUTPUT;
    @Argument(shortName = "REFFLAT", doc = "The refFlat gene model file")
    public File REFFLAT;
    @Argument(shortName = "T", doc = "The number of threads (default 20)")
    public int nThreads = 20;
    @Argument(shortName = "DELTA", doc = "Allowed base number difference between start/end of exons and read block position (default=10)")
    public int DELTA = 10;
    @Argument(shortName = "SOFT", doc = "Transcripts exons can be smaller than LongReadRecord exons (detection of specific alternative exons like flip/flop gria2 of Pkm1/Pkm2)")
    public boolean SOFT = false;

    public MoleculeConsensus() {
        log = Log.getInstance(MoleculeConsensus.class);
        pl = new htsjdk.samtools.util.ProgressLogger(log);
    }

    protected int doWork() {
        int nb = 0;
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(REFFLAT);

        UCSCRefFlatParser model = new UCSCRefFlatParser(REFFLAT);
        LongreadParser bam = new LongreadParser(INPUT);
        MoleculeDataset dataset = new MoleculeDataset(bam, model, DELTA, SOFT);

        dataset.callConsensus(OUTPUT, nThreads);

        return 0;
    }

    public static void main(String[] args) {
        System.exit(new MoleculeConsensus().instanceMain(args));
    }
}
