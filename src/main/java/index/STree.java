package index;


/** 
 * Suffix tree
 * 
*/
public class STree {
    public String T;
    public Node root;
    private Node actNode, insideNode;
    public int len = -1;
    private int actEdge = -1, actLength = 0, remainder = 0;
    private int end = -1;

    public STree(String strs) {
        this.T = strs + "#";
        this.buildTree();
    }

    public int edgelength(Node node) {
        return (node.leaf ? this.end : node.end) - node.start + 1;
    }

    private boolean walkdown(Node node) {
        int length = this.edgelength(node);
        if (this.actLength >= length) {
            this.actEdge += length;
            this.actLength -= length;
            this.actNode = node;
            return true;
        }
        return false;
    }

    private Node genNode(int start, int end, boolean leaf) {
        Node node = new Node(leaf);
        node.suffix_link = this.root;
        node.start = start;
        node.end = end;
        return node;
    }

    private void genTree(int pos) {
        this.end = pos;
        this.remainder += 1;
        this.insideNode = null;

        while (this.remainder > 0) {
            if (this.actLength == 0) {
                this.actEdge = pos;
            }
            if (this.actNode.children.get(this.T.charAt(this.actEdge)) == null) {
                this.actNode.children.put(this.T.charAt(this.actEdge), this.genNode(pos, -1, true));
                if (this.insideNode != null) {
                    this.insideNode.suffix_link = this.actNode;
                    this.insideNode = null;
                }
            }
            else {
                Node nextNode = this.actNode.children.get(this.T.charAt(this.actEdge));
                if (this.walkdown(nextNode)) {
                    continue;
                }
                if (this.T.charAt(nextNode.start + this.actLength) == this.T.charAt(pos)) {
                    if (this.insideNode != null && (this.actNode != this.root)) {
                        this.insideNode.suffix_link = this.actNode;
                        this.insideNode = null;
                    }

                    this.actLength += 1;
                    break;
                }
                int splitEnd = nextNode.start + this.actLength - 1;
                Node splitNode = this.genNode(nextNode.start, splitEnd, false);
                // this.actNode.children.remove(this.T.charAt(this.actEdge));
                this.actNode.children.put(this.T.charAt(this.actEdge), splitNode);
                splitNode.children.put(this.T.charAt(pos), this.genNode(pos, -1, true));
                nextNode.start += this.actLength;
                splitNode.children.put(this.T.charAt(nextNode.start), nextNode);

                if (this.insideNode != null) {
                    this.insideNode.suffix_link = splitNode; 
                }
                this.insideNode = splitNode;
            }
            this.remainder -= 1;
            if (this.actNode == this.root && this.remainder > 0) {
                this.actLength -= 1;
                this.actEdge = pos - this.remainder + 1;
            }
            else if (this.actNode != this.root) {
                this.actNode = this.actNode.suffix_link;
            }
        }
    }

    private void buildTree() {
        this.len = this.T.length();
        int rootEnd = -1;
        this.root = this.genNode(-1, rootEnd, false);
        this.actNode = this.root;
        for (int i = 0; i < this.len; ++i) {
            this.genTree(i);
        }
    }
}