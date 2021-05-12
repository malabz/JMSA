package main;

import java.io.IOException;
import io.fasta;
import msa.centerAlign;
import msa.treeAlign;

public class Main {
    private static String mode;
    private static String infile;
    private static String outfile;
    
    public static void main (String[] args) throws IOException {
        parse(args);
        print_args();
        // args_help();

        System.out.print("Read data START-->");
        fasta readfasta = new fasta();
        String[][] res = readfasta.readFasta(infile);
        String[] lables = res[0];
        String[] strs = res[1];
        System.out.println("DONE!");
        if (mode.equals("treealign")) {
            treeAlign talign = new treeAlign(strs);
            String[] strsTal = talign.getStrsAlign();
            readfasta.writeFasta(strsTal, lables, outfile);
        }
        else {
            centerAlign calign = new centerAlign(strs, lables, "suffix");
            String[] strsal = calign.getStrsAlign();
            readfasta.writeFasta(strsal, lables, outfile);
        }
    }

    private static void parse(String[] args) {
        // 比对方式选择 -m
        // 输入文件位置 -i
        // 输出文件位置 -o
        if (args.length == 0 || args.length > 6) { args_help(); }
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-m") && args.length > i + 1) {
                if (args[i+1].equalsIgnoreCase("centeralign") || args[i+1].equalsIgnoreCase("treealign")) {
                    mode = args[++i].toLowerCase();
                }
                else { 
                    System.out.println(args[i+1]);
                    args_help(); 
                }
            }
            else if (args[i].equals("-i") && args.length > i + 1) { infile = args[++i]; }
            else if (args[i].equals("-o") && args.length > i + 1) { outfile = args[++i]; }
            else { 
                args_help(); }
        }

        if (infile == null) args_help();
        if (mode == null) mode = "treealign";
        if (outfile == null) {
            if (mode.equals("treealign")) outfile = infile + ".TreeAligned";
            if (mode.equals("centeralign")) outfile = infile + ".CenterAligned";
        }    
    }

    private static void args_help()
    {
        System.out.println("\nusage: java -jar " + " [-m] mode [-i] path [-o] path" );
        System.out.println();
        System.out.println("  necessary arguments: ");
        System.out.println("    -i  Input file path (nucleotide sequences in fasta format)");
        System.out.println("    -m  two option (1.TreeAlign 2.CenterAlign)");
        System.out.println();
        System.out.println("  optional arguments: ");
        System.out.println("    -o  Output file path (default: inputfile.aligned)");
        System.out.println();
        System.exit(0);
    }

    private static void print_args()
    {
        System.out.println("\n**");
        System.out.println("** mode: " + mode);
        System.out.println("** infile: " + infile);
        System.out.println("** outfile: " + outfile);
        System.out.println("**\n");
        System.out.println();
    }
}
