package index;

import java.util.HashMap;

public class Node {
    
    public int start, end;
    public boolean leaf;
    public Node suffix_link;
    public HashMap<Character, Node> children = new HashMap<>();

    public Node (Boolean leaf) { this.leaf = leaf; }
}
