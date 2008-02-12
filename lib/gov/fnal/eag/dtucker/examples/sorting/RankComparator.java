package gov.fnal.eag.dtucker.examples.sorting;

//From http://java.sun.com/developer/JDCTechTips/2004/tt0716.html#1


import java.util.Comparator;

public class RankComparator implements Comparator {
	public int compare(Object o1, Object o2) {
		Rank r1 = (Rank) o1;
		Rank r2 = (Rank) o2;
		return r1.getRankInt() - r2.getRankInt();
	}
}

