/**
 * @file D2Q777.h
 * @short D2Q777 Stencil
 * @date 09.10.2018
 * @author dwilde3m, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef D2Q777MODEL_H_
#define D2Q777MODEL_H_

#include "Stencil.h"

#include "../utilities/BasicNames.h"

namespace natrium {

/** @short D2Q777 Model
 */
class D2Q777: public Stencil {

private:

	/** @short function to create the vector of directions
	 *  @return the vector of directions
	 */
	vector<numeric_vector> makeDirections(double scaling);

	/** @short function to create the vector of weights
	 *  @return the vector of weights
	 */
	vector<double> makeWeights();

	numeric_matrix makeMomentBasis(vector<numeric_vector> e);

protected:

	/// speed of sound
	const double m_speedOfSound;

	/// (speed of sound)^2
	const double m_speedOfSoundSquare;

	/// scaling of the stencil
	const double m_scaling;

public:

	/// D
	static constexpr size_t D=2;

	/// Q
	static constexpr size_t Q=777;

	/// constructor
	D2Q777(double scaling = 1.0);

	/// destructor
	virtual ~D2Q777();

	virtual double getSpeedOfSound() const {
		return m_speedOfSound;
	}
	virtual double getSpeedOfSoundSquare() const {
		return m_speedOfSoundSquare;
	}

	virtual size_t getIndexOfOppositeDirection(size_t index) const {
		switch (index){
            case 0 :
                return 776 ;
                break;
            case 1 :
                return 775 ;
                break;
            case 2 :
                return 774 ;
                break;
            case 3 :
                return 773 ;
                break;
            case 4 :
                return 772 ;
                break;
            case 5 :
                return 771 ;
                break;
            case 6 :
                return 770 ;
                break;
            case 7 :
                return 769 ;
                break;
            case 8 :
                return 768 ;
                break;
            case 9 :
                return 767 ;
                break;
            case 10 :
                return 766 ;
                break;
            case 11 :
                return 765 ;
                break;
            case 12 :
                return 764 ;
                break;
            case 13 :
                return 763 ;
                break;
            case 14 :
                return 762 ;
                break;
            case 15 :
                return 761 ;
                break;
            case 16 :
                return 760 ;
                break;
            case 17 :
                return 759 ;
                break;
            case 18 :
                return 758 ;
                break;
            case 19 :
                return 757 ;
                break;
            case 20 :
                return 756 ;
                break;
            case 21 :
                return 755 ;
                break;
            case 22 :
                return 754 ;
                break;
            case 23 :
                return 753 ;
                break;
            case 24 :
                return 752 ;
                break;
            case 25 :
                return 751 ;
                break;
            case 26 :
                return 750 ;
                break;
            case 27 :
                return 749 ;
                break;
            case 28 :
                return 748 ;
                break;
            case 29 :
                return 747 ;
                break;
            case 30 :
                return 746 ;
                break;
            case 31 :
                return 745 ;
                break;
            case 32 :
                return 744 ;
                break;
            case 33 :
                return 743 ;
                break;
            case 34 :
                return 742 ;
                break;
            case 35 :
                return 741 ;
                break;
            case 36 :
                return 740 ;
                break;
            case 37 :
                return 739 ;
                break;
            case 38 :
                return 738 ;
                break;
            case 39 :
                return 737 ;
                break;
            case 40 :
                return 736 ;
                break;
            case 41 :
                return 735 ;
                break;
            case 42 :
                return 734 ;
                break;
            case 43 :
                return 733 ;
                break;
            case 44 :
                return 732 ;
                break;
            case 45 :
                return 731 ;
                break;
            case 46 :
                return 730 ;
                break;
            case 47 :
                return 729 ;
                break;
            case 48 :
                return 728 ;
                break;
            case 49 :
                return 727 ;
                break;
            case 50 :
                return 726 ;
                break;
            case 51 :
                return 725 ;
                break;
            case 52 :
                return 724 ;
                break;
            case 53 :
                return 723 ;
                break;
            case 54 :
                return 722 ;
                break;
            case 55 :
                return 721 ;
                break;
            case 56 :
                return 720 ;
                break;
            case 57 :
                return 719 ;
                break;
            case 58 :
                return 718 ;
                break;
            case 59 :
                return 717 ;
                break;
            case 60 :
                return 716 ;
                break;
            case 61 :
                return 715 ;
                break;
            case 62 :
                return 714 ;
                break;
            case 63 :
                return 713 ;
                break;
            case 64 :
                return 712 ;
                break;
            case 65 :
                return 711 ;
                break;
            case 66 :
                return 710 ;
                break;
            case 67 :
                return 709 ;
                break;
            case 68 :
                return 708 ;
                break;
            case 69 :
                return 707 ;
                break;
            case 70 :
                return 706 ;
                break;
            case 71 :
                return 705 ;
                break;
            case 72 :
                return 704 ;
                break;
            case 73 :
                return 703 ;
                break;
            case 74 :
                return 702 ;
                break;
            case 75 :
                return 701 ;
                break;
            case 76 :
                return 700 ;
                break;
            case 77 :
                return 699 ;
                break;
            case 78 :
                return 698 ;
                break;
            case 79 :
                return 697 ;
                break;
            case 80 :
                return 696 ;
                break;
            case 81 :
                return 695 ;
                break;
            case 82 :
                return 694 ;
                break;
            case 83 :
                return 693 ;
                break;
            case 84 :
                return 692 ;
                break;
            case 85 :
                return 691 ;
                break;
            case 86 :
                return 690 ;
                break;
            case 87 :
                return 689 ;
                break;
            case 88 :
                return 688 ;
                break;
            case 89 :
                return 687 ;
                break;
            case 90 :
                return 686 ;
                break;
            case 91 :
                return 685 ;
                break;
            case 92 :
                return 684 ;
                break;
            case 93 :
                return 683 ;
                break;
            case 94 :
                return 682 ;
                break;
            case 95 :
                return 681 ;
                break;
            case 96 :
                return 680 ;
                break;
            case 97 :
                return 679 ;
                break;
            case 98 :
                return 678 ;
                break;
            case 99 :
                return 677 ;
                break;
            case 100 :
                return 676 ;
                break;
            case 101 :
                return 675 ;
                break;
            case 102 :
                return 674 ;
                break;
            case 103 :
                return 673 ;
                break;
            case 104 :
                return 672 ;
                break;
            case 105 :
                return 671 ;
                break;
            case 106 :
                return 670 ;
                break;
            case 107 :
                return 669 ;
                break;
            case 108 :
                return 668 ;
                break;
            case 109 :
                return 667 ;
                break;
            case 110 :
                return 666 ;
                break;
            case 111 :
                return 665 ;
                break;
            case 112 :
                return 664 ;
                break;
            case 113 :
                return 663 ;
                break;
            case 114 :
                return 662 ;
                break;
            case 115 :
                return 661 ;
                break;
            case 116 :
                return 660 ;
                break;
            case 117 :
                return 659 ;
                break;
            case 118 :
                return 658 ;
                break;
            case 119 :
                return 657 ;
                break;
            case 120 :
                return 656 ;
                break;
            case 121 :
                return 655 ;
                break;
            case 122 :
                return 654 ;
                break;
            case 123 :
                return 653 ;
                break;
            case 124 :
                return 652 ;
                break;
            case 125 :
                return 651 ;
                break;
            case 126 :
                return 650 ;
                break;
            case 127 :
                return 649 ;
                break;
            case 128 :
                return 648 ;
                break;
            case 129 :
                return 647 ;
                break;
            case 130 :
                return 646 ;
                break;
            case 131 :
                return 645 ;
                break;
            case 132 :
                return 644 ;
                break;
            case 133 :
                return 643 ;
                break;
            case 134 :
                return 642 ;
                break;
            case 135 :
                return 641 ;
                break;
            case 136 :
                return 640 ;
                break;
            case 137 :
                return 639 ;
                break;
            case 138 :
                return 638 ;
                break;
            case 139 :
                return 637 ;
                break;
            case 140 :
                return 636 ;
                break;
            case 141 :
                return 635 ;
                break;
            case 142 :
                return 634 ;
                break;
            case 143 :
                return 633 ;
                break;
            case 144 :
                return 632 ;
                break;
            case 145 :
                return 631 ;
                break;
            case 146 :
                return 630 ;
                break;
            case 147 :
                return 629 ;
                break;
            case 148 :
                return 628 ;
                break;
            case 149 :
                return 627 ;
                break;
            case 150 :
                return 626 ;
                break;
            case 151 :
                return 625 ;
                break;
            case 152 :
                return 624 ;
                break;
            case 153 :
                return 623 ;
                break;
            case 154 :
                return 622 ;
                break;
            case 155 :
                return 621 ;
                break;
            case 156 :
                return 620 ;
                break;
            case 157 :
                return 619 ;
                break;
            case 158 :
                return 618 ;
                break;
            case 159 :
                return 617 ;
                break;
            case 160 :
                return 616 ;
                break;
            case 161 :
                return 615 ;
                break;
            case 162 :
                return 614 ;
                break;
            case 163 :
                return 613 ;
                break;
            case 164 :
                return 612 ;
                break;
            case 165 :
                return 611 ;
                break;
            case 166 :
                return 610 ;
                break;
            case 167 :
                return 609 ;
                break;
            case 168 :
                return 608 ;
                break;
            case 169 :
                return 607 ;
                break;
            case 170 :
                return 606 ;
                break;
            case 171 :
                return 605 ;
                break;
            case 172 :
                return 604 ;
                break;
            case 173 :
                return 603 ;
                break;
            case 174 :
                return 602 ;
                break;
            case 175 :
                return 601 ;
                break;
            case 176 :
                return 600 ;
                break;
            case 177 :
                return 599 ;
                break;
            case 178 :
                return 598 ;
                break;
            case 179 :
                return 597 ;
                break;
            case 180 :
                return 596 ;
                break;
            case 181 :
                return 595 ;
                break;
            case 182 :
                return 594 ;
                break;
            case 183 :
                return 593 ;
                break;
            case 184 :
                return 592 ;
                break;
            case 185 :
                return 591 ;
                break;
            case 186 :
                return 590 ;
                break;
            case 187 :
                return 589 ;
                break;
            case 188 :
                return 588 ;
                break;
            case 189 :
                return 587 ;
                break;
            case 190 :
                return 586 ;
                break;
            case 191 :
                return 585 ;
                break;
            case 192 :
                return 584 ;
                break;
            case 193 :
                return 583 ;
                break;
            case 194 :
                return 582 ;
                break;
            case 195 :
                return 581 ;
                break;
            case 196 :
                return 580 ;
                break;
            case 197 :
                return 579 ;
                break;
            case 198 :
                return 578 ;
                break;
            case 199 :
                return 577 ;
                break;
            case 200 :
                return 576 ;
                break;
            case 201 :
                return 575 ;
                break;
            case 202 :
                return 574 ;
                break;
            case 203 :
                return 573 ;
                break;
            case 204 :
                return 572 ;
                break;
            case 205 :
                return 571 ;
                break;
            case 206 :
                return 570 ;
                break;
            case 207 :
                return 569 ;
                break;
            case 208 :
                return 568 ;
                break;
            case 209 :
                return 567 ;
                break;
            case 210 :
                return 566 ;
                break;
            case 211 :
                return 565 ;
                break;
            case 212 :
                return 564 ;
                break;
            case 213 :
                return 563 ;
                break;
            case 214 :
                return 562 ;
                break;
            case 215 :
                return 561 ;
                break;
            case 216 :
                return 560 ;
                break;
            case 217 :
                return 559 ;
                break;
            case 218 :
                return 558 ;
                break;
            case 219 :
                return 557 ;
                break;
            case 220 :
                return 556 ;
                break;
            case 221 :
                return 555 ;
                break;
            case 222 :
                return 554 ;
                break;
            case 223 :
                return 553 ;
                break;
            case 224 :
                return 552 ;
                break;
            case 225 :
                return 551 ;
                break;
            case 226 :
                return 550 ;
                break;
            case 227 :
                return 549 ;
                break;
            case 228 :
                return 548 ;
                break;
            case 229 :
                return 547 ;
                break;
            case 230 :
                return 546 ;
                break;
            case 231 :
                return 545 ;
                break;
            case 232 :
                return 544 ;
                break;
            case 233 :
                return 543 ;
                break;
            case 234 :
                return 542 ;
                break;
            case 235 :
                return 541 ;
                break;
            case 236 :
                return 540 ;
                break;
            case 237 :
                return 539 ;
                break;
            case 238 :
                return 538 ;
                break;
            case 239 :
                return 537 ;
                break;
            case 240 :
                return 536 ;
                break;
            case 241 :
                return 535 ;
                break;
            case 242 :
                return 534 ;
                break;
            case 243 :
                return 533 ;
                break;
            case 244 :
                return 532 ;
                break;
            case 245 :
                return 531 ;
                break;
            case 246 :
                return 530 ;
                break;
            case 247 :
                return 529 ;
                break;
            case 248 :
                return 528 ;
                break;
            case 249 :
                return 527 ;
                break;
            case 250 :
                return 526 ;
                break;
            case 251 :
                return 525 ;
                break;
            case 252 :
                return 524 ;
                break;
            case 253 :
                return 523 ;
                break;
            case 254 :
                return 522 ;
                break;
            case 255 :
                return 521 ;
                break;
            case 256 :
                return 520 ;
                break;
            case 257 :
                return 519 ;
                break;
            case 258 :
                return 518 ;
                break;
            case 259 :
                return 517 ;
                break;
            case 260 :
                return 516 ;
                break;
            case 261 :
                return 515 ;
                break;
            case 262 :
                return 514 ;
                break;
            case 263 :
                return 513 ;
                break;
            case 264 :
                return 512 ;
                break;
            case 265 :
                return 511 ;
                break;
            case 266 :
                return 510 ;
                break;
            case 267 :
                return 509 ;
                break;
            case 268 :
                return 508 ;
                break;
            case 269 :
                return 507 ;
                break;
            case 270 :
                return 506 ;
                break;
            case 271 :
                return 505 ;
                break;
            case 272 :
                return 504 ;
                break;
            case 273 :
                return 503 ;
                break;
            case 274 :
                return 502 ;
                break;
            case 275 :
                return 501 ;
                break;
            case 276 :
                return 500 ;
                break;
            case 277 :
                return 499 ;
                break;
            case 278 :
                return 498 ;
                break;
            case 279 :
                return 497 ;
                break;
            case 280 :
                return 496 ;
                break;
            case 281 :
                return 495 ;
                break;
            case 282 :
                return 494 ;
                break;
            case 283 :
                return 493 ;
                break;
            case 284 :
                return 492 ;
                break;
            case 285 :
                return 491 ;
                break;
            case 286 :
                return 490 ;
                break;
            case 287 :
                return 489 ;
                break;
            case 288 :
                return 488 ;
                break;
            case 289 :
                return 487 ;
                break;
            case 290 :
                return 486 ;
                break;
            case 291 :
                return 485 ;
                break;
            case 292 :
                return 484 ;
                break;
            case 293 :
                return 483 ;
                break;
            case 294 :
                return 482 ;
                break;
            case 295 :
                return 481 ;
                break;
            case 296 :
                return 480 ;
                break;
            case 297 :
                return 479 ;
                break;
            case 298 :
                return 478 ;
                break;
            case 299 :
                return 477 ;
                break;
            case 300 :
                return 476 ;
                break;
            case 301 :
                return 475 ;
                break;
            case 302 :
                return 474 ;
                break;
            case 303 :
                return 473 ;
                break;
            case 304 :
                return 472 ;
                break;

            case 305 :
                return 471 ;
                break;
            case 306 :
                return 470 ;
                break;
            case 307 :
                return 469 ;
                break;
            case 308 :
                return 468 ;
                break;
            case 309 :
                return 467 ;
                break;
            case 310 :
                return 466 ;
                break;
            case 311 :
                return 465 ;
                break;
            case 312 :
                return 464 ;
                break;
            case 313 :
                return 463 ;
                break;
            case 314 :
                return 462 ;
                break;
            case 315 :
                return 461 ;
                break;
            case 316 :
                return 460 ;
                break;
            case 317 :
                return 459 ;
                break;
            case 318 :
                return 458 ;
                break;
            case 319 :
                return 457 ;
                break;
            case 320 :
                return 456 ;
                break;
            case 321 :
                return 455 ;
                break;
            case 322 :
                return 454 ;
                break;
            case 323 :
                return 453 ;
                break;
            case 324 :
                return 452 ;
                break;
            case 325 :
                return 451 ;
                break;
            case 326 :
                return 450 ;
                break;
            case 327 :
                return 449 ;
                break;
            case 328 :
                return 448 ;
                break;
            case 329 :
                return 447 ;
                break;
            case 330 :
                return 446 ;
                break;
            case 331 :
                return 445 ;
                break;
            case 332 :
                return 444 ;
                break;
            case 333 :
                return 443 ;
                break;
            case 334 :
                return 442 ;
                break;
            case 335 :
                return 441 ;
                break;
            case 336 :
                return 440 ;
                break;
            case 337 :
                return 439 ;
                break;
            case 338 :
                return 438 ;
                break;
            case 339 :
                return 437 ;
                break;
            case 340 :
                return 436 ;
                break;
            case 341 :
                return 435 ;
                break;
            case 342 :
                return 434 ;
                break;
            case 343 :
                return 433 ;
                break;
            case 344 :
                return 432 ;
                break;
            case 345 :
                return 431 ;
                break;
            case 346 :
                return 430 ;
                break;
            case 347 :
                return 429 ;
                break;
            case 348 :
                return 428 ;
                break;
            case 349 :
                return 427 ;
                break;
            case 350 :
                return 426 ;
                break;
            case 351 :
                return 425 ;
                break;
            case 352 :
                return 424 ;
                break;
            case 353 :
                return 423 ;
                break;
            case 354 :
                return 422 ;
                break;
            case 355 :
                return 421 ;
                break;
            case 356 :
                return 420 ;
                break;
            case 357 :
                return 419 ;
                break;
            case 358 :
                return 418 ;
                break;
            case 359 :
                return 417 ;
                break;
            case 360 :
                return 416 ;
                break;
            case 361 :
                return 415 ;
                break;
            case 362 :
                return 414 ;
                break;
            case 363 :
                return 413 ;
                break;
            case 364 :
                return 412 ;
                break;
            case 365 :
                return 411 ;
                break;
            case 366 :
                return 410 ;
                break;
            case 367 :
                return 409 ;
                break;
            case 368 :
                return 408 ;
                break;
            case 369 :
                return 407 ;
                break;
            case 370 :
                return 406 ;
                break;
            case 371 :
                return 405 ;
                break;
            case 372 :
                return 404 ;
                break;
            case 373 :
                return 403 ;
                break;
            case 374 :
                return 402 ;
                break;
            case 375 :
                return 401 ;
                break;
            case 376 :
                return 400 ;
                break;
            case 377 :
                return 399 ;
                break;
            case 378 :
                return 398 ;
                break;
            case 379 :
                return 397 ;
                break;
            case 380 :
                return 396 ;
                break;
            case 381 :
                return 395 ;
                break;
            case 382 :
                return 394 ;
                break;
            case 383 :
                return 393 ;
                break;
            case 384 :
                return 392 ;
                break;
            case 385 :
                return 391 ;
                break;
            case 386 :
                return 390 ;
                break;
            case 387 :
                return 389 ;
                break;
            case 388 :
                return 388 ;
                break;
            case 389 :
                return 387 ;
                break;
            case 390 :
                return 386 ;
                break;
            case 391 :
                return 385 ;
                break;
            case 392 :
                return 384 ;
                break;
            case 393 :
                return 383 ;
                break;
            case 394 :
                return 382 ;
                break;
            case 395 :
                return 381 ;
                break;
            case 396 :
                return 380 ;
                break;
            case 397 :
                return 379 ;
                break;
            case 398 :
                return 378 ;
                break;
            case 399 :
                return 377 ;
                break;
            case 400 :
                return 376 ;
                break;
            case 401 :
                return 375 ;
                break;
            case 402 :
                return 374 ;
                break;
            case 403 :
                return 373 ;
                break;
            case 404 :
                return 372 ;
                break;
            case 405 :
                return 371 ;
                break;
            case 406 :
                return 370 ;
                break;
            case 407 :
                return 369 ;
                break;
            case 408 :
                return 368 ;
                break;
            case 409 :
                return 367 ;
                break;
            case 410 :
                return 366 ;
                break;
            case 411 :
                return 365 ;
                break;
            case 412 :
                return 364 ;
                break;
            case 413 :
                return 363 ;
                break;
            case 414 :
                return 362 ;
                break;
            case 415 :
                return 361 ;
                break;
            case 416 :
                return 360 ;
                break;
            case 417 :
                return 359 ;
                break;
            case 418 :
                return 358 ;
                break;
            case 419 :
                return 357 ;
                break;
            case 420 :
                return 356 ;
                break;
            case 421 :
                return 355 ;
                break;
            case 422 :
                return 354 ;
                break;
            case 423 :
                return 353 ;
                break;
            case 424 :
                return 352 ;
                break;
            case 425 :
                return 351 ;
                break;
            case 426 :
                return 350 ;
                break;
            case 427 :
                return 349 ;
                break;
            case 428 :
                return 348 ;
                break;
            case 429 :
                return 347 ;
                break;
            case 430 :
                return 346 ;
                break;
            case 431 :
                return 345 ;
                break;
            case 432 :
                return 344 ;
                break;
            case 433 :
                return 343 ;
                break;
            case 434 :
                return 342 ;
                break;
            case 435 :
                return 341 ;
                break;
            case 436 :
                return 340 ;
                break;
            case 437 :
                return 339 ;
                break;
            case 438 :
                return 338 ;
                break;
            case 439 :
                return 337 ;
                break;
            case 440 :
                return 336 ;
                break;
            case 441 :
                return 335 ;
                break;
            case 442 :
                return 334 ;
                break;
            case 443 :
                return 333 ;
                break;
            case 444 :
                return 332 ;
                break;
            case 445 :
                return 331 ;
                break;
            case 446 :
                return 330 ;
                break;
            case 447 :
                return 329 ;
                break;
            case 448 :
                return 328 ;
                break;
            case 449 :
                return 327 ;
                break;
            case 450 :
                return 326 ;
                break;
            case 451 :
                return 325 ;
                break;
            case 452 :
                return 324 ;
                break;
            case 453 :
                return 323 ;
                break;
            case 454 :
                return 322 ;
                break;
            case 455 :
                return 321 ;
                break;
            case 456 :
                return 320 ;
                break;
            case 457 :
                return 319 ;
                break;
            case 458 :
                return 318 ;
                break;
            case 459 :
                return 317 ;
                break;
            case 460 :
                return 316 ;
                break;
            case 461 :
                return 315 ;
                break;
            case 462 :
                return 314 ;
                break;
            case 463 :
                return 313 ;
                break;
            case 464 :
                return 312 ;
                break;
            case 465 :
                return 311 ;
                break;
            case 466 :
                return 310 ;
                break;
            case 467 :
                return 309 ;
                break;
            case 468 :
                return 308 ;
                break;
            case 469 :
                return 307 ;
                break;
            case 470 :
                return 306 ;
                break;
            case 471 :
                return 305 ;
                break;
            case 472 :
                return 304 ;
                break;
            case 473 :
                return 303 ;
                break;
            case 474 :
                return 302 ;
                break;
            case 475 :
                return 301 ;
                break;
            case 476 :
                return 300 ;
                break;
            case 477 :
                return 299 ;
                break;
            case 478 :
                return 298 ;
                break;
            case 479 :
                return 297 ;
                break;
            case 480 :
                return 296 ;
                break;
            case 481 :
                return 295 ;
                break;
            case 482 :
                return 294 ;
                break;
            case 483 :
                return 293 ;
                break;
            case 484 :
                return 292 ;
                break;
            case 485 :
                return 291 ;
                break;
            case 486 :
                return 290 ;
                break;
            case 487 :
                return 289 ;
                break;
            case 488 :
                return 288 ;
                break;
            case 489 :
                return 287 ;
                break;
            case 490 :
                return 286 ;
                break;
            case 491 :
                return 285 ;
                break;
            case 492 :
                return 284 ;
                break;
            case 493 :
                return 283 ;
                break;
            case 494 :
                return 282 ;
                break;
            case 495 :
                return 281 ;
                break;
            case 496 :
                return 280 ;
                break;
            case 497 :
                return 279 ;
                break;
            case 498 :
                return 278 ;
                break;
            case 499 :
                return 277 ;
                break;
            case 500 :
                return 276 ;
                break;
            case 501 :
                return 275 ;
                break;
            case 502 :
                return 274 ;
                break;
            case 503 :
                return 273 ;
                break;
            case 504 :
                return 272 ;
                break;
            case 505 :
                return 271 ;
                break;
            case 506 :
                return 270 ;
                break;
            case 507 :
                return 269 ;
                break;
            case 508 :
                return 268 ;
                break;
            case 509 :
                return 267 ;
                break;
            case 510 :
                return 266 ;
                break;
            case 511 :
                return 265 ;
                break;
            case 512 :
                return 264 ;
                break;
            case 513 :
                return 263 ;
                break;
            case 514 :
                return 262 ;
                break;
            case 515 :
                return 261 ;
                break;
            case 516 :
                return 260 ;
                break;
            case 517 :
                return 259 ;
                break;
            case 518 :
                return 258 ;
                break;
            case 519 :
                return 257 ;
                break;
            case 520 :
                return 256 ;
                break;
            case 521 :
                return 255 ;
                break;
            case 522 :
                return 254 ;
                break;
            case 523 :
                return 253 ;
                break;
            case 524 :
                return 252 ;
                break;
            case 525 :
                return 251 ;
                break;
            case 526 :
                return 250 ;
                break;
            case 527 :
                return 249 ;
                break;
            case 528 :
                return 248 ;
                break;
            case 529 :
                return 247 ;
                break;
            case 530 :
                return 246 ;
                break;
            case 531 :
                return 245 ;
                break;
            case 532 :
                return 244 ;
                break;
            case 533 :
                return 243 ;
                break;
            case 534 :
                return 242 ;
                break;
            case 535 :
                return 241 ;
                break;
            case 536 :
                return 240 ;
                break;
            case 537 :
                return 239 ;
                break;
            case 538 :
                return 238 ;
                break;
            case 539 :
                return 237 ;
                break;
            case 540 :
                return 236 ;
                break;
            case 541 :
                return 235 ;
                break;
            case 542 :
                return 234 ;
                break;
            case 543 :
                return 233 ;
                break;
            case 544 :
                return 232 ;
                break;
            case 545 :
                return 231 ;
                break;
            case 546 :
                return 230 ;
                break;
            case 547 :
                return 229 ;
                break;
            case 548 :
                return 228 ;
                break;
            case 549 :
                return 227 ;
                break;
            case 550 :
                return 226 ;
                break;
            case 551 :
                return 225 ;
                break;
            case 552 :
                return 224 ;
                break;
            case 553 :
                return 223 ;
                break;
            case 554 :
                return 222 ;
                break;
            case 555 :
                return 221 ;
                break;
            case 556 :
                return 220 ;
                break;
            case 557 :
                return 219 ;
                break;
            case 558 :
                return 218 ;
                break;
            case 559 :
                return 217 ;
                break;
            case 560 :
                return 216 ;
                break;
            case 561 :
                return 215 ;
                break;
            case 562 :
                return 214 ;
                break;
            case 563 :
                return 213 ;
                break;
            case 564 :
                return 212 ;
                break;
            case 565 :
                return 211 ;
                break;
            case 566 :
                return 210 ;
                break;
            case 567 :
                return 209 ;
                break;
            case 568 :
                return 208 ;
                break;
            case 569 :
                return 207 ;
                break;
            case 570 :
                return 206 ;
                break;
            case 571 :
                return 205 ;
                break;
            case 572 :
                return 204 ;
                break;
            case 573 :
                return 203 ;
                break;
            case 574 :
                return 202 ;
                break;
            case 575 :
                return 201 ;
                break;
            case 576 :
                return 200 ;
                break;
            case 577 :
                return 199 ;
                break;
            case 578 :
                return 198 ;
                break;
            case 579 :
                return 197 ;
                break;
            case 580 :
                return 196 ;
                break;
            case 581 :
                return 195 ;
                break;
            case 582 :
                return 194 ;
                break;
            case 583 :
                return 193 ;
                break;
            case 584 :
                return 192 ;
                break;
            case 585 :
                return 191 ;
                break;
            case 586 :
                return 190 ;
                break;
            case 587 :
                return 189 ;
                break;
            case 588 :
                return 188 ;
                break;
            case 589 :
                return 187 ;
                break;
            case 590 :
                return 186 ;
                break;
            case 591 :
                return 185 ;
                break;
            case 592 :
                return 184 ;
                break;
            case 593 :
                return 183 ;
                break;
            case 594 :
                return 182 ;
                break;
            case 595 :
                return 181 ;
                break;
            case 596 :
                return 180 ;
                break;
            case 597 :
                return 179 ;
                break;
            case 598 :
                return 178 ;
                break;
            case 599 :
                return 177 ;
                break;
            case 600 :
                return 176 ;
                break;
            case 601 :
                return 175 ;
                break;
            case 602 :
                return 174 ;
                break;
            case 603 :
                return 173 ;
                break;
            case 604 :
                return 172 ;
                break;
            case 605 :
                return 171 ;
                break;
            case 606 :
                return 170 ;
                break;
            case 607 :
                return 169 ;
                break;
            case 608 :
                return 168 ;
                break;

            case 609 :
                return 167 ;
                break;
            case 610 :
                return 166 ;
                break;
            case 611 :
                return 165 ;
                break;
            case 612 :
                return 164 ;
                break;
            case 613 :
                return 163 ;
                break;
            case 614 :
                return 162 ;
                break;
            case 615 :
                return 161 ;
                break;
            case 616 :
                return 160 ;
                break;
            case 617 :
                return 159 ;
                break;
            case 618 :
                return 158 ;
                break;
            case 619 :
                return 157 ;
                break;
            case 620 :
                return 156 ;
                break;
            case 621 :
                return 155 ;
                break;
            case 622 :
                return 154 ;
                break;
            case 623 :
                return 153 ;
                break;
            case 624 :
                return 152 ;
                break;
            case 625 :
                return 151 ;
                break;
            case 626 :
                return 150 ;
                break;
            case 627 :
                return 149 ;
                break;
            case 628 :
                return 148 ;
                break;
            case 629 :
                return 147 ;
                break;
            case 630 :
                return 146 ;
                break;
            case 631 :
                return 145 ;
                break;
            case 632 :
                return 144 ;
                break;
            case 633 :
                return 143 ;
                break;
            case 634 :
                return 142 ;
                break;
            case 635 :
                return 141 ;
                break;
            case 636 :
                return 140 ;
                break;
            case 637 :
                return 139 ;
                break;
            case 638 :
                return 138 ;
                break;
            case 639 :
                return 137 ;
                break;
            case 640 :
                return 136 ;
                break;
            case 641 :
                return 135 ;
                break;
            case 642 :
                return 134 ;
                break;
            case 643 :
                return 133 ;
                break;
            case 644 :
                return 132 ;
                break;
            case 645 :
                return 131 ;
                break;
            case 646 :
                return 130 ;
                break;
            case 647 :
                return 129 ;
                break;
            case 648 :
                return 128 ;
                break;
            case 649 :
                return 127 ;
                break;
            case 650 :
                return 126 ;
                break;
            case 651 :
                return 125 ;
                break;
            case 652 :
                return 124 ;
                break;
            case 653 :
                return 123 ;
                break;
            case 654 :
                return 122 ;
                break;
            case 655 :
                return 121 ;
                break;
            case 656 :
                return 120 ;
                break;
            case 657 :
                return 119 ;
                break;
            case 658 :
                return 118 ;
                break;
            case 659 :
                return 117 ;
                break;
            case 660 :
                return 116 ;
                break;
            case 661 :
                return 115 ;
                break;
            case 662 :
                return 114 ;
                break;
            case 663 :
                return 113 ;
                break;
            case 664 :
                return 112 ;
                break;
            case 665 :
                return 111 ;
                break;
            case 666 :
                return 110 ;
                break;
            case 667 :
                return 109 ;
                break;
            case 668 :
                return 108 ;
                break;
            case 669 :
                return 107 ;
                break;
            case 670 :
                return 106 ;
                break;
            case 671 :
                return 105 ;
                break;
            case 672 :
                return 104 ;
                break;
            case 673 :
                return 103 ;
                break;
            case 674 :
                return 102 ;
                break;
            case 675 :
                return 101 ;
                break;
            case 676 :
                return 100 ;
                break;
            case 677 :
                return 99 ;
                break;
            case 678 :
                return 98 ;
                break;
            case 679 :
                return 97 ;
                break;
            case 680 :
                return 96 ;
                break;
            case 681 :
                return 95 ;
                break;
            case 682 :
                return 94 ;
                break;
            case 683 :
                return 93 ;
                break;
            case 684 :
                return 92 ;
                break;
            case 685 :
                return 91 ;
                break;
            case 686 :
                return 90 ;
                break;
            case 687 :
                return 89 ;
                break;
            case 688 :
                return 88 ;
                break;
            case 689 :
                return 87 ;
                break;
            case 690 :
                return 86 ;
                break;
            case 691 :
                return 85 ;
                break;
            case 692 :
                return 84 ;
                break;
            case 693 :
                return 83 ;
                break;
            case 694 :
                return 82 ;
                break;
            case 695 :
                return 81 ;
                break;
            case 696 :
                return 80 ;
                break;
            case 697 :
                return 79 ;
                break;
            case 698 :
                return 78 ;
                break;
            case 699 :
                return 77 ;
                break;
            case 700 :
                return 76 ;
                break;
            case 701 :
                return 75 ;
                break;
            case 702 :
                return 74 ;
                break;
            case 703 :
                return 73 ;
                break;
            case 704 :
                return 72 ;
                break;
            case 705 :
                return 71 ;
                break;
            case 706 :
                return 70 ;
                break;
            case 707 :
                return 69 ;
                break;
            case 708 :
                return 68 ;
                break;
            case 709 :
                return 67 ;
                break;
            case 710 :
                return 66 ;
                break;
            case 711 :
                return 65 ;
                break;
            case 712 :
                return 64 ;
                break;
            case 713 :
                return 63 ;
                break;
            case 714 :
                return 62 ;
                break;
            case 715 :
                return 61 ;
                break;
            case 716 :
                return 60 ;
                break;
            case 717 :
                return 59 ;
                break;
            case 718 :
                return 58 ;
                break;
            case 719 :
                return 57 ;
                break;
            case 720 :
                return 56 ;
                break;
            case 721 :
                return 55 ;
                break;
            case 722 :
                return 54 ;
                break;
            case 723 :
                return 53 ;
                break;
            case 724 :
                return 52 ;
                break;
            case 725 :
                return 51 ;
                break;
            case 726 :
                return 50 ;
                break;
            case 727 :
                return 49 ;
                break;
            case 728 :
                return 48 ;
                break;
            case 729 :
                return 47 ;
                break;
            case 730 :
                return 46 ;
                break;
            case 731 :
                return 45 ;
                break;
            case 732 :
                return 44 ;
                break;
            case 733 :
                return 43 ;
                break;
            case 734 :
                return 42 ;
                break;
            case 735 :
                return 41 ;
                break;
            case 736 :
                return 40 ;
                break;
            case 737 :
                return 39 ;
                break;
            case 738 :
                return 38 ;
                break;
            case 739 :
                return 37 ;
                break;
            case 740 :
                return 36 ;
                break;
            case 741 :
                return 35 ;
                break;
            case 742 :
                return 34 ;
                break;
            case 743 :
                return 33 ;
                break;
            case 744 :
                return 32 ;
                break;
            case 745 :
                return 31 ;
                break;
            case 746 :
                return 30 ;
                break;
            case 747 :
                return 29 ;
                break;
            case 748 :
                return 28 ;
                break;
            case 749 :
                return 27 ;
                break;
            case 750 :
                return 26 ;
                break;
            case 751 :
                return 25 ;
                break;
            case 752 :
                return 24 ;
                break;
            case 753 :
                return 23 ;
                break;
            case 754 :
                return 22 ;
                break;
            case 755 :
                return 21 ;
                break;
            case 756 :
                return 20 ;
                break;
            case 757 :
                return 19 ;
                break;
            case 758 :
                return 18 ;
                break;
            case 759 :
                return 17 ;
                break;
            case 760 :
                return 16 ;
                break;
            case 761 :
                return 15 ;
                break;
            case 762 :
                return 14 ;
                break;
            case 763 :
                return 13 ;
                break;
            case 764 :
                return 12 ;
                break;
            case 765 :
                return 11 ;
                break;
            case 766 :
                return 10 ;
                break;
            case 767 :
                return 9 ;
                break;
            case 768 :
                return 8 ;
                break;
            case 769 :
                return 7 ;
                break;
            case 770 :
                return 6 ;
                break;
            case 771 :
                return 5 ;
                break;
            case 772 :
                return 4 ;
                break;
            case 773 :
                return 3 ;
                break;
            case 774 :
                return 2 ;
                break;
            case 775 :
                return 1 ;
                break;
            case 776 :
                return 0 ;
                break;

		}
	}

	virtual double getMaxParticleVelocityMagnitude() const {
        return sqrt(2)*m_scaling;
	}

	virtual double getScaling() const {
		return m_scaling;
	}

};

} /* namespace natrium */
#endif /* D2Q777MODEL_H_ */
