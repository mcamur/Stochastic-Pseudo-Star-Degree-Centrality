
public class Edge {
	private int node1;
	private int node2;
	public int getNode1() {
		return node1;
	}
	public void setNode1(int node1) {
		this.node1 = node1;
	}
	public int getNode2() {
		return node2;
	}
	public void setNode2(int node2) {
		this.node2 = node2;
	}
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + node1;
		result = prime * result + node2;
		return result;
	}
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		Edge other = (Edge) obj;
		if (node1 != other.node1)
			return false;
		if (node2 != other.node2)
			return false;
		return true;
	}
	public Edge(int node1, int node2) {
		super();
		this.node1 = node1;
		this.node2 = node2;
	}
	@Override
	public String toString() {
		return "(" + node1 + "," + node2 + ")";
	}

	
	
}
