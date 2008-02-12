package gov.fnal.eag.dtucker.examples.sorting;


//From http://java.sun.com/developer/JDCTechTips/2004/tt0716.html#1

public class UseComparableRank extends ShuffleAndSort {
	
	void createCard(int i) {
		miniDeck.add(new ComparableRank(i));
	}
	
	public static void main(String[] args) {
		new UseComparableRank().exerciseDeck();
	}
}

