package msa;

/**
 * cenTree 还需要好好挖掘一下潜力，比如说到底是该
 * 不相似的序列是不是不适合用suffix
 * 还是适合的，但是如何选择呢
 * 第一遍先用centeralign
 * 第二遍用treealign
 */
public class CenTree {
    private String[] straligned;
    private final String Treemode;

    public CenTree ( String[] strs, String treemode) {
        this.Treemode = treemode;
        Align(strs);
    }

    public String[] getStrsAlign() {
        return this.straligned;
    }

    private void Align(String[] strs) {
        centerAlign cAlign = new centerAlign(strs, 1);
        straligned = cAlign.getStrsAlign();
        treeAlign tAlign = new treeAlign(strs, straligned, Treemode);
        straligned = tAlign.getStrsAlign();
    }
}