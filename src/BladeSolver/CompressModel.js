//Web worker for compressing/decompressing stringified models.
/*
	http://pieroxy.net/blog/pages/lz-string			
	var string = "This is my compression test."; 
	alert("Size of sample is: " + string.length); 
	var compressed = LZString.compress(string); 
	alert("Size of compressed sample is: " + compressed.length); 
	string = LZString.decompress(compressed); 
	alert("Sample is: " + string);			
*/

self.onmessage = function(e){
	self.importScripts('../js/lz-string.min.js');

	if(e.data.plaintext){
		var compressedText = LZString.compressToBase64(e.data.plaintext);
		var message = ("plaintext length="+e.data.plaintext.length+" compressed length="+compressedText.length);
		self.postMessage({"compressedText":compressedText, "message":message});		
	}
	else if(e.data.compressedText){
		var plaintext = LZString.decompressFromBase64(e.data.compressedText);
		var message = ("plaintext length="+plaintext.length+" compressed length="+e.data.compressedText.length);		
		self.postMessage({"plaintext":plaintext, "message":message});		
	}
	
	close();
};
