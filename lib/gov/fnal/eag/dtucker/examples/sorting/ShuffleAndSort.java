package gov.fnal.eag.dtucker.examples.sorting;

//From http://java.sun.com/developer/JDCTechTips/2004/tt0716.html#1

import java.util.List;
import java.util.ArrayList;
import java.util.Collections;

public class ShuffleAndSort {
	List miniDeck = new ArrayList(6);
	
	void initializeDeck() {
		for (int i = 0; i < 6; i++) {
			createCard(i);
		}
	}
	
	void printDeck(String message) {
		System.out.println(message);
		for (int i = 0; i < 6; i++) {
			System.out.println("card " + i +
					" = " + miniDeck.get(i));
		}
		System.out.println("============");
	}
	
	void createCard(int i) {
		miniDeck.add(new Integer(i));
	}
	
	void sort() {
		Collections.sort(miniDeck);
	}
	
	//Sort in reverse order...
	//void sort(){
    //    Collections.sort(miniDeck,
	//                       Collections.reverseOrder());
	//}  

	
	void exerciseDeck() {
		initializeDeck();
		Collections.shuffle(miniDeck);
		printDeck("Deck Shuffled:");
		sort();
		printDeck("Deck Sorted:");
	}
	
	
	public static void main(String[] args) {
		new ShuffleAndSort().exerciseDeck();
	}
}
