package gov.fnal.eag.dtucker.examples.sorting;

//From http://java.sun.com/developer/JDCTechTips/2004/tt0716.html#1

public class Rank {
	
	private int rank;
	
	public Rank(int rank) {
		this.rank = rank;
	}
	
	public int getRankInt(){
		return rank;
	}
	
	public String toString() {
		return "" + rank;
	}
}

