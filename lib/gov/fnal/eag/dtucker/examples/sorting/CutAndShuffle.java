package gov.fnal.eag.dtucker.examples.sorting;

//From http://java.sun.com/developer/JDCTechTips/2004/tt0716.html#1

import java.util.List;
import java.util.ArrayList;
import java.util.Collections;

public class CutAndShuffle {   
	
	private static List miniDeck = new ArrayList(6);
	
	private static void initializeDeck() {
		for (int i = 0; i < 6; i++) {
			miniDeck.add(new Integer(i));
		}
	}
	
	private static void printDeck(String message) {
		System.out.println(message + "\n");
		for (int i = 0; i < 6; i++) {
			System.out.println("card " + i +
					" = " + miniDeck.get(i));
		}
		System.out.println("============");
	}
	
	public static void main(String[] args) {
		initializeDeck();
		printDeck("Initialized Deck:");
		Collections.rotate(miniDeck, 2);
		printDeck("Deck rotated by 2:");
		Collections.rotate(miniDeck, -2);
		printDeck("Deck rotated back by 2:");
		Collections.shuffle(miniDeck);
		printDeck("Deck Shuffled:");
	}
}

