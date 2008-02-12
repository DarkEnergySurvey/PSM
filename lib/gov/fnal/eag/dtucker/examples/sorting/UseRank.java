package gov.fnal.eag.dtucker.examples.sorting;

//From http://java.sun.com/developer/JDCTechTips/2004/tt0716.html#1

//The main method of this class will fail, since Rank not does implement Comparable.
//When you compile and run UseRank you get a ClassCastException. 
//There are two ways to solve this problem:
//    * You create a Comparable version of Rank.
//    * You store the rank as an int and implement compareTo() "by hand."

public class UseRank extends ShuffleAndSort {
	
	void createCard(int i) {
		miniDeck.add(new Rank(i));
	}
	
	public static void main(String[] args) {
		new UseRank().exerciseDeck();
	}
}

