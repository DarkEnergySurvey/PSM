package gov.fnal.eag.dtucker.examples.sorting;

//From http://java.sun.com/developer/JDCTechTips/2004/tt0716.html#1

public class ComparableRank2 implements Comparable {
	
	private int rank;
	
	public ComparableRank2(int rank) {
		this.rank = rank;
	}
	
	public String toString() {
		return "" + rank;
	}
	
	public int compareTo(Object o) {
		ComparableRank2 cr = (ComparableRank2) o;
		return (rank - cr.rank);
	}
}

