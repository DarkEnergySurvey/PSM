package gov.fnal.eag.dtucker.examples.sorting;

//From http://java.sun.com/developer/JDCTechTips/2004/tt0716.html#1

public class ComparableRank implements Comparable {
	
	private Integer rank;
	
	public ComparableRank(int rank) {
		this.rank = new Integer(rank);
	}
	
	public String toString() {
		return "" + rank;
	}
	
	public int compareTo(Object o) {
		ComparableRank cr = (ComparableRank) o;
		return (rank.compareTo(cr.rank));
	}
}
